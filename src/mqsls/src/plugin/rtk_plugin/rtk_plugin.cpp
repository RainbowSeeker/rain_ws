#include <rclcpp/rclcpp.hpp>
#include <boost/asio.hpp>
#include <mqsls/msg/follower_send.hpp>
#include <iostream>
#include <string>
#include <array>

#include "rtk_plugin.hpp"
#include "mqsls/perf_counter.hpp"
#include "mqsls/SlidingWindowFilter.hpp"

#define RTK_UNICORE     0
#define RTK_SEPTENTRIO  1
#define RTK_MANUFACTURER_SELECT    RTK_SEPTENTRIO

#if RTK_MANUFACTURER_SELECT == RTK_UNICORE
#include "unicore.hpp"
#elif RTK_MANUFACTURER_SELECT == RTK_SEPTENTRIO
#include "septentrio.hpp"
#else
#error "RTK_MANUFACTURER_SELECT is not defined"
#endif

using namespace boost::asio;

class RTKPlugin : public rclcpp::Node
{
public:
    RTKPlugin() : Node("rtk_plugin")
    {
        // parameters check
        if (_amc_id < 0 || _amc_id > 3) {
            throw std::invalid_argument("Invalid amc_id");
        }

        // configure serial port
        _serial.open(_port_name);
        _serial.set_option(serial_port::baud_rate(_baudrate));
        _serial.set_option(serial_port::flow_control(serial_port::flow_control::none));
        _serial.set_option(serial_port::parity(serial_port::parity::none));
        _serial.set_option(serial_port::stop_bits(serial_port::stop_bits::one));
        _serial.set_option(serial_port::character_size(8));

        // publisher
        _follower_send_pub = this->create_publisher<mqsls::msg::FollowerSend>("follower_send" + std::to_string(_amc_id), 10);

        RCLCPP_INFO(this->get_logger(), "RTK plugin started for AMC %d", _amc_id);
    }

    ~RTKPlugin() {}

    void start()
    {
        // select component
#if RTK_MANUFACTURER_SELECT == RTK_UNICORE
        _rtk_component = std::make_shared<unicore::DualAntennaComponent>(this->shared_from_this());
#elif RTK_MANUFACTURER_SELECT == RTK_SEPTENTRIO
        if (_amc_id == 0) {
            _rtk_component = std::make_shared<septentrio::BaseStationComponent>(this->shared_from_this());
        } else {
            _rtk_component = std::make_shared<septentrio::RoverComponent>(this->shared_from_this(), [this](const std::vector<uint8_t> &bytes) {
                _serial.async_write_some(buffer(bytes), [=](const boost::system::error_code &error, size_t bytes_transferred) {
                    if (error)
                    {
                        RCLCPP_ERROR(this->get_logger(), "Error writing to serial port: %s", error.message().c_str());
                    }
                });
            });
        }
#endif

        do_async_read();

        std::thread([this]() {
            _io.run();
        }).detach();
    }

private:
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
        // handle the message
        auto result = RefTranslation();
        int ret = _rtk_component->handle_recv_message(raw, result);
        if (ret > 0) {
            // publish the result
            auto msg = mqsls::msg::FollowerSend();
            msg.position_uav[0] = result.delta_xyz[0];
            msg.position_uav[1] = result.delta_xyz[1];
            msg.position_uav[2] = result.delta_xyz[2];
            msg.timestamp = this->get_clock()->now().nanoseconds() / 1e3; // [us]
            _follower_send_pub->publish(msg);

            _perf_counter.tick();
            if (_perf_counter.count() % 100 == 0) {
                _perf_counter.print();
            }
        }
    }
    
    // parameters
    const int _amc_id = this->declare_parameter("amc_id", 1); // 0 : base station; 1-3 : rover
    const std::string _port_name = this->declare_parameter<std::string>("port_name", "/dev/ttyUSB0");
    const int _baudrate = this->declare_parameter<int>("baudrate", 115200);

    // publisher
    rclcpp::Publisher<mqsls::msg::FollowerSend>::SharedPtr _follower_send_pub;

    // details
    RTKComponent::SharedPtr _rtk_component;

    // filter
    SlidingWindowFilter<RefTranslation> _filter {10};
    PerfCounter _perf_counter {PerfCounterType::INTERVAL, "RTKPlugin"};

    // serial port
    io_context _io;
    serial_port _serial {_io};
    std::array<char, 10240> _recv_buf;
};

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);

    auto node = std::make_shared<RTKPlugin>();

    node->start();
    
    rclcpp::spin(node);
    
    rclcpp::shutdown();

    return 0;
}