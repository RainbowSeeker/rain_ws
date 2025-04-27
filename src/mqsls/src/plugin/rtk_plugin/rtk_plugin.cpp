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

#ifndef DEG2RAD
#define DEG2RAD(x) ((x) * 0.01745329251994329575)
#endif

#ifndef RAD2DEG
#define RAD2DEG(x) ((x) * 57.29577951308232087721)
#endif

using namespace boost::asio;

class RTKPlugin : public rclcpp::Node
{
public:
    RTKPlugin() : Node("rtk_plugin")
    {
        // configure serial port
        _serial.open(_port_name);
        _serial.set_option(serial_port::baud_rate(_baudrate));
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
                        _follower_send_msg.timestamp = absolute_time();
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
        _parser.handle_message(raw, [this]() {
            auto msg_id = _parser.get_message_id();
            switch (msg_id)
            {
                case septentrio::msg::MSG_ID_POSITION_CARTESIAN:
                {
                    auto body = _parser.get_position_cartesian();
                    _follower_send_msg.position_load[0] = body.X;
                    _follower_send_msg.position_load[1] = body.Y;
                    _follower_send_msg.position_load[2] = body.Z;
                    _follower_send_msg.position_uav[0] = body.Base2RoverX;
                    _follower_send_msg.position_uav[1] = body.Base2RoverY;
                    _follower_send_msg.position_uav[2] = body.Base2RoverZ;
                    _follower_send_msg.timestamp = absolute_time();
                    RCLCPP_INFO(this->get_logger(), "POSITION_CARTESIAN: %f, %f, %f", _follower_send_msg.position_uav[0], _follower_send_msg.position_uav[1], _follower_send_msg.position_uav[2]);
                    break;
                }
                case septentrio::msg::MSG_ID_ATTITUDE_EULER:
                {
                    auto body = _parser.get_attitude_euler();
                    double heading_rad = DEG2RAD(body.Heading);
                    if (heading_rad > M_PI) {
                        heading_rad -= 2. * M_PI;
                    }
                    double pitch_rad = DEG2RAD(body.Pitch);

                    // auto filtered = _filter.update({body.Heading, heading_rad, pitch_rad});

                    // // Update uav position && load position
                    // _follower_send_msg.position_uav[0] = 0;
                    // _follower_send_msg.position_uav[1] = 0;
                    // _follower_send_msg.position_uav[2] = 0;
                    // _follower_send_msg.position_load[0] = _follower_send_msg.position_uav[0] + cosf(heading_rad) * cosf(pitch_rad) * body.Base2RoverX;
                    // _follower_send_msg.position_load[1] = _follower_send_msg.position_uav[1] + sinf(heading_rad) * cosf(pitch_rad) * body.Base2RoverY;
                    // _follower_send_msg.position_load[2] = _follower_send_msg.position_uav[2] - sinf(pitch_rad) * body.Base2RoverZ;
                    // _follower_send_msg.timestamp = absolute_time();
                    // _follower_send_pub->publish(_follower_send_msg);

                    RCLCPP_INFO(this->get_logger(), "ATTITUDE_EULER: heading: %.2f, pitch: %.2f, roll: %.2f", body.Heading, body.Pitch, body.Roll);
                    RCLCPP_INFO(this->get_logger(), "ATTITUDE_EULER: heading_dot: %.2f, pitch_dot: %.2f, roll_dot: %.2f", body.HeadingDot, body.PitchDot, body.RollDot);

                    _perf_counter.tick();
                    if (_perf_counter.count() % 100 == 0) {
                        _perf_counter.print();
                    }
                    break;
                }
                default:
                    RCLCPP_ERROR(this->get_logger(), "Unknown message ID: %d", msg_id);
                    break;
            }
        });
#endif
    }

    inline uint64_t absolute_time() {
        return this->get_clock()->now().nanoseconds() / 1e3; // [us]
    }
    
    // parameters
    const int _amc_id = this->declare_parameter("amc_id", 1);

    // decode
#if RTK_MANUFACTURER_SELECT == RTK_SEPTENTRIO
    septentrio::msg::Parser _parser {};
#endif

    // publisher
    rclcpp::Publisher<mqsls::msg::FollowerSend>::SharedPtr _follower_send_pub;
    mqsls::msg::FollowerSend _follower_send_msg {};

    // status
    std::atomic<bool> _is_running = false;
    PerfCounter _perf_counter {PerfCounterType::INTERVAL, "RTKPlugin"};

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

    // serial
    const std::string _port_name = this->declare_parameter<std::string>("port_name", "/dev/ttyUSB0");
    const int _baudrate = this->declare_parameter<int>("baudrate", 115200);
    io_context _io;
    serial_port _serial {_io};
    std::array<char, 1024> _recv_buf;
};

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<RTKPlugin>());
    
    rclcpp::shutdown();

    return 0;
}