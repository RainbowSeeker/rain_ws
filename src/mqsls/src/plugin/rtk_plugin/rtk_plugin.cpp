#include <rclcpp/rclcpp.hpp>
#include <boost/asio.hpp>
#include <mqsls/msg/follower_send.hpp>
#include <iostream>
#include <string>
#include <array>

#include "mqsls/perf_counter.hpp"
#include "unicore.hpp"

using namespace boost::asio;

class RTKPlugin : public rclcpp::Node
{
public:
    RTKPlugin(const std::string &port_name = "/dev/ttyUSB0", const int &baudrate = 115200)
     : Node("rtk_plugin"), _serial(_io, port_name)
    {
        // configure serial port
        _serial.set_option(serial_port::baud_rate(baudrate));
        _serial.set_option(serial_port::flow_control(serial_port::flow_control::none));
        _serial.set_option(serial_port::parity(serial_port::parity::none));
        _serial.set_option(serial_port::stop_bits(serial_port::stop_bits::one));
        _serial.set_option(serial_port::character_size(8));

        // ros2
        _follower_send_pub = this->create_publisher<mqsls::msg::FollowerSend>("follower_send" + std::to_string(_amc_id), 10);

        // start
        start();

        RCLCPP_INFO(this->get_logger(), "RTK plugin started for AMC %d", _amc_id);
    }

    ~RTKPlugin() {}

private:
    void start()
    {
        if (_is_running) {
            RCLCPP_INFO(this->get_logger(), "Already running");
            return;
        }
        
        do_async_read();

        std::thread([this]() {
            _io.run();
        }).detach();

        _is_running = true;
    }

    void do_async_read()
    {
        _serial.async_read_some(buffer(_recv_buf), [this](const boost::system::error_code &error, size_t bytes_transferred) {
            if (!error)
            {
                std::string recv_data(_recv_buf.data(), bytes_transferred);
                parse_package(recv_data);

                do_async_read();
            }
            else
            {
                RCLCPP_ERROR(this->get_logger(), "Error reading serial port: %s", error.message().c_str());
            }
        });
    }

    void parse_package(const std::string &raw)
    {
        static unicore::msg::body_binary msg = {};
        for (const auto &ch : raw)
        {
            if (unicore::parse(ch, msg) == 0) {
                switch (msg.header.msg_id) {
                    case unicore::msg::MSG_ID_BESTNAVXYZ:
                    {
                        auto body = reinterpret_cast<unicore::msg::bestnavxyz *>(&msg.payload);
                        _follower_send_msg.velocity_uav[0] = body->vel[0];
                        _follower_send_msg.velocity_uav[1] = -body->vel[1];
                        _follower_send_msg.velocity_uav[2] = -body->vel[2];
                        // RCLCPP_INFO(this->get_logger(), "BESTNAVXYZ: %f, %f, %f", _follower_send_msg.velocity_uav[0], _follower_send_msg.velocity_uav[1], _follower_send_msg.velocity_uav[2]);
                        break;
                    }
                    case unicore::msg::MSG_ID_BESTNAVXYZH:
                    {
                        auto body = reinterpret_cast<unicore::msg::bestnavxyz *>(&msg.payload);
                        _follower_send_msg.velocity_load[0] = body->vel[0];
                        _follower_send_msg.velocity_load[1] = -body->vel[1];
                        _follower_send_msg.velocity_load[2] = -body->vel[2];
                        // RCLCPP_INFO(this->get_logger(), "BESTNAVXYZH: %f, %f, %f", _follower_send_msg.velocity_load[0], _follower_send_msg.velocity_load[1], _follower_send_msg.velocity_load[2]);
                        break;
                    }
                    case unicore::msg::MSG_ID_UNIHEADING:
                    {
                        auto body = reinterpret_cast<unicore::msg::uniheading *>(&msg.payload);
                        double heading_rad = body->heading * M_PI / 180.0;
                        if (heading_rad > M_PI) {
                            heading_rad -= 2. * M_PI;
                        }
                        double pitch_rad = body->pitch * M_PI / 180.0;
        
                        // Update uav position && load position
                        _follower_send_msg.position_uav[0] = 0;
                        _follower_send_msg.position_uav[1] = 0;
                        _follower_send_msg.position_uav[2] = 0;
                        _follower_send_msg.position_load[0] = _follower_send_msg.position_uav[0] + cosf(heading_rad) * cosf(pitch_rad) * body->baseline;
                        _follower_send_msg.position_load[1] = _follower_send_msg.position_uav[1] + sinf(heading_rad) * cosf(pitch_rad) * body->baseline;
                        _follower_send_msg.position_load[2] = _follower_send_msg.position_uav[2] - sinf(pitch_rad) * body->baseline;
                        _follower_send_msg.timestamp = absolute_time();
                        _follower_send_pub->publish(_follower_send_msg);
        
                        // RCLCPP_INFO(this->get_logger(), "UNIHEADING: len: %.3f, heading: %.2f, pitch: %.2f", (double)body->baseline, (double)body->heading, (double)body->pitch);
                        // RCLCPP_INFO(this->get_logger(), "UNIHEADING: delta_pos: %.3f, %.3f, %.3f", _follower_send_msg.position_load[0] - _follower_send_msg.position_uav[0], _follower_send_msg.position_load[1] - _follower_send_msg.position_uav[1], _follower_send_msg.position_load[2] - _follower_send_msg.position_uav[2]);

                        _perf_counter.tick();
                        if (_perf_counter.count() % 100 == 0) {
                            _perf_counter.print();
                        }
                        break;
                    }
                    default:
                        RCLCPP_ERROR(this->get_logger(), "Unknown message ID: %d", msg.header.msg_id);
                        break;
                }
            }
        }
    }

    inline uint64_t absolute_time() {
        return this->get_clock()->now().nanoseconds() / 1e3; // [us]
    }
    
    // parameters
    const int _amc_id = this->declare_parameter("amc_id", 1);

    // publisher
    rclcpp::Publisher<mqsls::msg::FollowerSend>::SharedPtr _follower_send_pub;
    mqsls::msg::FollowerSend _follower_send_msg {};

    // status
    std::atomic<bool> _is_running = false;
    PerfCounter _perf_counter {PerfCounterType::INTERVAL, "RTKPlugin"};

    // serial
    io_context _io;
    serial_port _serial;
    std::array<char, 1024> _recv_buf;
};

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<RTKPlugin>());
    
    rclcpp::shutdown();

    return 0;
}