#include <boost/asio.hpp>
#include <rclcpp/rclcpp.hpp>
#include <Eigen/Eigen>
#include <iostream>

#include <px4_msgs/msg/vehicle_odometry.hpp>
#include <mqsls/msg/follower_send.hpp>

using namespace boost::asio;

#define SAMPLE_FREQ     100 // Hz

struct MocapData
{
    int id; // 1, 2, 3 for uav; 4 for payload
    uint64_t timestamp;
    float position[3];  // Motion capture data
    float q[4];
};

class MocapPlugin : public rclcpp::Node
{
public:
    MocapPlugin() : rclcpp::Node("mocap_plugin"),
        _udp_endpoint(ip::udp::v4(), 14555), _udp_socket(_io, _udp_endpoint)
    {
        RCLCPP_INFO(this->get_logger(), "MocapPlugin node started");

        const std::string topic_ns = "/px4_";

        // Publisher
        // Only 1, 2, 3 are used for PX4 && Leader
        for (int i = 0; i < 3; i++)
        {
            _odom_pub[i] = this->create_publisher<px4_msgs::msg::VehicleOdometry>(
                topic_ns + std::to_string(i + 1) + "/fmu/in/vehicle_visual_odometry", 10);

            _follower_pub[i] = this->create_publisher<mqsls::msg::FollowerSend>(
                "/follower_send" + std::to_string(i + 1), 10);
        }

        start_connection();
    }

    ~MocapPlugin()
    {
        _io.stop();
    }

private:
    void start_connection()
    {
        do_udp_recv();

        std::thread([this]() {
            _io.run();
        }).detach();
    }

    void do_udp_recv()
    {
        _udp_socket.async_receive_from(
            boost::asio::buffer(_recv_buf), _remote_endpoint,
            [this](const boost::system::error_code &error, std::size_t bytes_transferred) {
                if (!error)
                {
                    handle_recv_message(_recv_buf.data(), bytes_transferred);
                    do_udp_recv();
                }
                else
                {
                    std::cerr << "Error receiving data: " << error.message() << std::endl;
                }
            });
    }


    void handle_recv_message(const char *buf, size_t len)
    {
        if (len != sizeof(struct MocapData))
        {
            RCLCPP_ERROR(this->get_logger(), "Invalid data size: %lu", len);
            return;
        }

        auto data = reinterpret_cast<const struct MocapData *>(buf);

        if (data->id < 1 || data->id > 4)
        {
            RCLCPP_ERROR(this->get_logger(), "Invalid ID: %d", data->id);
            return;
        }

        // To PX4
        if (data->id <= 3)
        {
            // To PX4
            px4_msgs::msg::VehicleOdometry odom {};
            odom.timestamp = data->timestamp;
            odom.position[0] = data->position[0];
            odom.position[1] = data->position[1];
            odom.position[2] = data->position[2];
            odom.q[0] = data->q[0];
            odom.q[1] = data->q[1];
            odom.q[2] = data->q[2];
            odom.q[3] = data->q[3];

            _odom_pub[data->id - 1]->publish(odom);
        }

        // To Leader
        if (data->id <= 3)
        {
            for (int i = 0; i < 3; i++)
            {
                _uav_velocity[data->id - 1][i] = (data->position[i] - _uav_position[data->id - 1][i]) * SAMPLE_FREQ;
                _uav_position[data->id - 1][i] = data->position[i];
            }

            mqsls::msg::FollowerSend send {};
            send.timestamp = data->timestamp;
            for (int i = 0; i < 3; i++)
            {
                send.position_load[i] = _load_positon[i];
                send.velocity_load[i] = _load_velocity[i];
                send.position_uav[i] = _uav_position[data->id - 1][i];
                send.velocity_uav[i] = _uav_velocity[data->id - 1][i];
            }
            _follower_pub[data->id - 1]->publish(send);
        }
        else
        {
            for (int i = 0; i < 3; i++)
            {
                _load_velocity[i] = (data->position[i] - _load_positon[i]) * SAMPLE_FREQ;
                _load_positon[i] = data->position[i];
            }
        }


        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, 
            "ID: %d, Pos: (%.2f, %.2f, %.2f)",
            data->id, data->position[0], data->position[1], data->position[2]);    
    }

    // asio
    io_context _io;
    ip::udp::endpoint _udp_endpoint;
    ip::udp::endpoint _remote_endpoint;
    ip::udp::socket _udp_socket;
    std::array<char, 1024> _recv_buf;

    // Publisher
    rclcpp::Publisher<px4_msgs::msg::VehicleOdometry>::SharedPtr _odom_pub[3];
    rclcpp::Publisher<mqsls::msg::FollowerSend>::SharedPtr _follower_pub[3];

    // state
    Eigen::Vector3f _load_positon;
    Eigen::Vector3f _load_velocity;
    Eigen::Vector3f _uav_position[3];
    Eigen::Vector3f _uav_velocity[3];
};


int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<MocapPlugin>());
    rclcpp::shutdown();
    return 0;
}