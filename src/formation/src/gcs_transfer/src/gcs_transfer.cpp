#include <iostream>
#include <thread>
#include <vector>
#include <chrono>
#include <boost/asio.hpp>
#include <rclcpp/rclcpp.hpp>

#include <mavlink.h>
#include <formation/msg/uav_command.hpp>
#include <formation/msg/uav_status.hpp>

using namespace std::chrono_literals;
using namespace formation::msg;
using boost::asio::ip::tcp;

class MavlinkSession : public std::enable_shared_from_this<MavlinkSession>
{
public:
    explicit MavlinkSession(std::shared_ptr<tcp::socket> socket, int uav_num)
        : _socket(socket), 
        _endpoint_name(socket->remote_endpoint().address().to_string()), 
        _uav_num(uav_num)
    {
    }

    ~MavlinkSession()
    {
        std::cout << "Connection to " << _endpoint_name << " closed" << std::endl;
    }

    void start()
    {
        start_node();

        start_read();
    }

    std::shared_ptr<rclcpp::Node> get_node() const
    {
        return _node;
    }

private:
    static std::string current_time()
    {
        auto now = std::chrono::system_clock::now();
        auto now_c = std::chrono::system_clock::to_time_t(now);
        auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
        std::stringstream ss;
        ss << std::put_time(std::localtime(&now_c), "%H_%M_%S");
        return ss.str() + "_" + std::to_string(now_ms.count());
    }

    void start_node()
    {
        _node = std::make_shared<rclcpp::Node>("mavlink_session_" + current_time());

        const char *topic_ns = "/px4_";

        _uav_status_sub.resize(_uav_num);
        _uav_command_pub.resize(_uav_num);
        for (int i = 0; i < _uav_num; i++)
        {
            _uav_status_sub[i] = _node->create_subscription<UavStatus>(
                topic_ns + std::to_string(i + 1) + "/fmu/out/uav_status", 10,
                [this, i](const UavStatus::SharedPtr sp) {
                    mavlink_message_t msg;
                    mavlink_uav_status_t status;
                    status.timestamp = sp->timestamp;
                    status.stage = sp->stage;
                    static_assert(sizeof(status.sta_msg) == sizeof(sp->sta_msg), "Size of sta_msg does not match");
                    std::memcpy(status.sta_msg, sp->sta_msg.data(), sp->sta_msg.size());
                    mavlink_msg_uav_status_encode(i + 1, 0, &msg, &status);
                    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
                    int size = mavlink_msg_to_send_buffer(buf, &msg);
                    if (size > 0)
                    {
                        boost::system::error_code error;
                        _socket->write_some(boost::asio::buffer(buf, size), error);
                    }
                });

            _uav_command_pub[i] = _node->create_publisher<UavCommand>(
                topic_ns + std::to_string(i + 1) + "/fmu/in/uav_command", 10);
        }    
    }

    void start_read()
    {
        auto self(shared_from_this());
        _socket->async_read_some(boost::asio::buffer(_read_buf), [this, self](const boost::system::error_code &ec, std::size_t len) {
            if (!ec)
            {
                for (std::size_t i = 0; i < len; i++)
                {
                    if (mavlink_parse_char(MAVLINK_COMM_0, _read_buf[i], &_mavlink_msg, &_mavlink_status))
                    {
                        handle_mavlink_message(_mavlink_msg);
                    }
                }
                start_read();
            }
        });
    }

    void handle_mavlink_message(const mavlink_message_t &msg)
    {
        if (msg.sysid < 1 || msg.sysid > _uav_num)
        {
            std::cerr << "Invalid sysid: " << msg.sysid << std::endl;
            return;
        }

        switch (msg.msgid)
        {
        case MAVLINK_MSG_ID_UAV_COMMAND:
            {
                mavlink_uav_command_t cmd;
                mavlink_msg_uav_command_decode(&msg, &cmd);

                auto uav_command = UavCommand();
                uav_command.timestamp = cmd.timestamp;
                uav_command.param1 = cmd.param1;
                uav_command.param2 = cmd.param2;
                uav_command.param3 = cmd.param3;
                uav_command.command = cmd.command;
                _uav_command_pub[msg.sysid - 1]->publish(uav_command);
                std::cout << "Received UAV_COMMAND message from UAV [" << (int)msg.sysid << "]" << " cmd: " << cmd.command << std::endl;
                break;
            }
        default:
            std::cerr << "Unknown message ID: " << msg.msgid << std::endl;
            break;
        }
    }

    std::shared_ptr<tcp::socket> _socket;
    const std::string _endpoint_name;
    std::array<char, 1024> _read_buf;
    std::shared_ptr<rclcpp::Node> _node;

    // Mavlink parser
    mavlink_message_t _mavlink_msg;
    mavlink_status_t _mavlink_status;
    
    // Part ros2
    int _uav_num;

    // Subscriber
    std::vector<rclcpp::Subscription<UavStatus>::SharedPtr> _uav_status_sub;

    // Publisher
    std::vector<rclcpp::Publisher<UavCommand>::SharedPtr> _uav_command_pub;
};


class GcsTransManager : public rclcpp::Node
{
public:
    explicit GcsTransManager(int uav_num) :
        Node("GcsTransManager"), _acceptor(_io_context), _uav_num(uav_num)
    {
        start_tcp_server();
    }

    ~GcsTransManager()
    {
        _io_context.stop();
        _io_thread.join();
    }

private:
    // Part boost::asio
    void start_tcp_server()
    {
        try
        {
            tcp::endpoint endpoint(tcp::v4(), 14440);
            _acceptor.open(endpoint.protocol());
            _acceptor.set_option(tcp::acceptor::reuse_address(true));
            _acceptor.bind(endpoint);
            _acceptor.listen();

            // Start the io_context in a separate thread
            _io_thread = std::thread([this]() {
                RCLCPP_INFO(get_logger(), "Starting TCP server on port 14440");
                _io_context.run();
            });

            // Accept incoming connections
            handle_accept();
        }
        catch (const std::exception &e)
        {
            RCLCPP_ERROR(get_logger(), "Error starting TCP server: %s", e.what());
        }
    }

    void handle_accept()
    {
        auto socket = std::make_shared<tcp::socket>(_io_context);
        _acceptor.async_accept(*socket, [this, socket](const boost::system::error_code &ec) {
            if (!ec)
            {
                RCLCPP_INFO(get_logger(), "New connection from %s", socket->remote_endpoint().address().to_string().c_str());

                std::thread([this, socket]() {
                    auto session = std::make_shared<MavlinkSession>(socket, _uav_num);
                    session->start();

                    rclcpp::Rate rate(100);
                    while (rclcpp::ok() && session.use_count() > 1) // Wait until the session is destroyed
                    {
                        rclcpp::spin_some(session->get_node());
                        rate.sleep();
                    }
                }).detach();
            }
            else
            {
                RCLCPP_ERROR(get_logger(), "Error accepting connection: %s", ec.message().c_str());
            }
            handle_accept();
        });
    }
    boost::asio::io_context _io_context;
    tcp::acceptor _acceptor;
    std::thread _io_thread;

    int _uav_num;
};

int main(int argc, char * argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: GcsTransManager <uav_num>" << std::endl;
        return 1;
    }

    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<GcsTransManager>(std::stoi(argv[1])));
    rclcpp::shutdown();
    return 0;
}
