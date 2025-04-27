#include <rclcpp/rclcpp.hpp>
#include <boost/asio.hpp>
#include <mqsls/msg/follower_send.hpp>
#include <iostream>
#include <string>
#include <array>

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

        RCLCPP_INFO(this->get_logger(), "RTK plugin started for AMC %d", _amc_id);
    }

    ~RTKPlugin() {}

    void start()
    {
        // select component
#if RTK_MANUFACTURER_SELECT == RTK_SEPTENTRIO
        if (_amc_id == 0) {
            _sep_component = std::make_shared<septentrio::BaseStationComponent>(this->shared_from_this());
        } else {
            _sep_component = std::make_shared<septentrio::RoverComponent>(this->shared_from_this(), [this](const std::vector<uint8_t> &bytes) {
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
#if RTK_MANUFACTURER_SELECT == RTK_UNICORE
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
                        RCLCPP_DEBUG(this->get_logger(), "BESTNAVXYZ: %f, %f, %f", _follower_send_msg.velocity_uav[0], _follower_send_msg.velocity_uav[1], _follower_send_msg.velocity_uav[2]);
                        break;
                    }
                    case unicore::msg::MSG_ID_BESTNAVXYZH:
                    {
                        auto body = reinterpret_cast<unicore::msg::bestnavxyz *>(&msg.payload);
                        _follower_send_msg.velocity_load[0] = body->vel[0];
                        _follower_send_msg.velocity_load[1] = -body->vel[1];
                        _follower_send_msg.velocity_load[2] = -body->vel[2];
                        RCLCPP_DEBUG(this->get_logger(), "BESTNAVXYZH: %f, %f, %f", _follower_send_msg.velocity_load[0], _follower_send_msg.velocity_load[1], _follower_send_msg.velocity_load[2]);
                        break;
                    }
                    case unicore::msg::MSG_ID_UNIHEADING:
                    {
                        auto body = reinterpret_cast<unicore::msg::uniheading *>(&msg.payload);
                        double heading_rad = DEG2RAD(body->heading);
                        if (heading_rad > M_PI) {
                            heading_rad -= 2. * M_PI;
                        }
                        double pitch_rad = DEG2RAD(body->pitch);

                        auto filtered = _filter.update({body->baseline, heading_rad, pitch_rad});
        
                        // Update uav position && load position
                        _follower_send_msg.position_uav[0] = 0;
                        _follower_send_msg.position_uav[1] = 0;
                        _follower_send_msg.position_uav[2] = 0;
                        _follower_send_msg.position_load[0] = _follower_send_msg.position_uav[0] + cosf(heading_rad) * cosf(pitch_rad) * body->baseline;
                        _follower_send_msg.position_load[1] = _follower_send_msg.position_uav[1] + sinf(heading_rad) * cosf(pitch_rad) * body->baseline;
                        _follower_send_msg.position_load[2] = _follower_send_msg.position_uav[2] - sinf(pitch_rad) * body->baseline;
                        _follower_send_msg.timestamp = this->get_clock()->now().nanoseconds() / 1e3; // [us]
                        _follower_send_pub->publish(_follower_send_msg);

                        RCLCPP_INFO(this->get_logger(), "UNIHEADING: len: %.3f, heading: %.2f, pitch: %.2f", filtered.baseline, RAD2DEG(filtered.heading), RAD2DEG(filtered.pitch));
                        RCLCPP_INFO(this->get_logger(), "UNIHEADING: delta_pos: %.3f, %.3f, %.3f", filtered.delta_x(), filtered.delta_y(), filtered.delta_z());

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
#elif RTK_MANUFACTURER_SELECT == RTK_SEPTENTRIO
        _sep_component->handle_recv_message(raw);
#endif
    }
    
    // parameters
    const int _amc_id = this->declare_parameter("amc_id", 1); // 0 : base station; 1-3 : rover
    const std::string _port_name = this->declare_parameter<std::string>("port_name", "/dev/ttyUSB0");
    const int _baudrate = this->declare_parameter<int>("baudrate", 115200);

#if RTK_MANUFACTURER_SELECT == RTK_SEPTENTRIO
    septentrio::Component::SharedPtr _sep_component;
#endif

    // filter
    struct RefTranslation
    {
        double baseline;
        double heading;
        double pitch;

        double delta_x() const
        {
            return baseline * cosf(heading) * cosf(pitch);
        }

        double delta_y() const
        {
            return baseline * sinf(heading) * cosf(pitch);
        }

        double delta_z() const
        {
            return -baseline * sinf(pitch);
        }

        const RefTranslation operator+(const RefTranslation &other) const
        {
            return {baseline + other.baseline, heading + other.heading, pitch + other.pitch};
        }
        const RefTranslation operator-(const RefTranslation &other) const
        {
            return {baseline - other.baseline, heading - other.heading, pitch - other.pitch};
        }
        const RefTranslation operator/(const double &other) const
        {
            return {baseline / other, heading / other, pitch / other};
        }
    };
    SlidingWindowFilter<RefTranslation> _filter {10};

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