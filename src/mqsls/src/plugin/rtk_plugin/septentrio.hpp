#pragma once
#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/byte_multi_array.hpp>

#include <vector>
#include <cstdint>
#include <string>
#include <functional>

#include "crc16.h"

namespace septentrio {
namespace msg
{
#pragma pack(push, 1)
    struct header_binary
    {
        uint8_t     sync[2]; 	// 0x24, 0x40
        uint16_t    crc16;	    // CRC16
        uint16_t    id;	        // message id: bits 0-12: block number; bits 13-15: block revision number
        uint16_t    len;	    // message length
    };

    struct PositionCartesian
    {
        uint32_t    TOW;        // GPS time of week [s]
        uint16_t    WNc;        // GPS week number [w]
        uint8_t     Type : 6;   // 0: no fix; 1: 2D fix; 2: 3D fix; 3: RTK float; 4: RTK fixed; 5: DR; 6: PVT
        uint8_t     AutoBase : 1; // 0: no auto base; 1: auto base
        uint8_t     Flag2D : 1;     // 2D/3D flag
        uint8_t     Error;      // PVT error code
        double      X;          // X coordinate in coordinate frame specified by Datum [m]
        double      Y;          // Y coordinate in coordinate frame specified by Datum [m]
        double      Z;          // Z coordinate in coordinate frame specified by Datum [m]
        double      Base2RoverX; // X baseline component (from base to rover) [m]
        double      Base2RoverY; // Y baseline component (from base to rover) [m]
        double      Base2RoverZ; // Z baseline component (from base to rover) [m]
        float       Cov_xx;      // Variance of the x estimate
        float       Cov_yy;      // Variance of the y estimate
        float       Cov_zz;      // Variance of the z estimate
        float       Cov_xy;      // Covariance of x and y
        float       Cov_xz;      // Covariance of x and z
        float       Cov_yz;      // Covariance of y and z
        uint16_t    PDOP;        // [0.01] If 0, PDOP not available, otherwise divide by 100 to obtain PDOP.
        uint16_t    HDOP;        // [0.01] If 0, HDOP not available, otherwise divide by 100 to obtain HDOP.
        uint16_t    VDOP;        // [0.01] If 0, VDOP not available, otherwise divide by 100 to obtain VDOP.
        uint8_t     Misc;
        uint8_t     Reserved[13];
    };

    struct AttitudeEuler
    {
        uint32_t    TOW;        // GPS time of week [s]
        uint16_t    WNc;        // GPS week number [w]
        uint8_t     NrSV;       // NumSV
        uint8_t     Error;      // PVT error code
        uint16_t    Mode;       // 0: no attitude
        uint16_t    Reserved;
        float       Heading;    // Heading [deg]
        float       Pitch;      // Pitch [deg]
        float       Roll;       // Roll [deg]
        float       PitchDot;   // Pitch rate [deg/s]
        float       RollDot;    // Roll rate [deg/s]
        float       HeadingDot; // Heading rate [deg/s]
    };

    static_assert(sizeof(header_binary) == 8, "Header size mismatch");
    static_assert(sizeof(PositionCartesian) == 108 - sizeof(header_binary), "PositionCartesian size mismatch");
    static_assert(sizeof(AttitudeEuler) == 44 - sizeof(header_binary), "AttitudeEuler size mismatch");
    

#pragma pack(pop)

    // Message ID
    enum MessageID : uint16_t
    {
        MSG_ID_POSITION_CARTESIAN = 4044,
        MSG_ID_ATTITUDE_EULER = 5938,
    };

    static const size_t MAX_MESSAGE_SIZE = 256;
    class Parser
    {
    public:
        Parser() {
            _buffer.reserve(MAX_MESSAGE_SIZE);
        }
        ~Parser() = default;

        /*
         * @brief 
         * @param bytes
         * @return 0: no message; >0: number of messages
         */
        void handle_message(const std::string &bytes, std::function<void()> callback) {
            for (const auto &ch : bytes) {
                handle_one_byte(ch, callback);
            }
        }

        MessageID get_message_id() const {
            if (_buffer.size() < sizeof(header_binary)) {
                throw std::runtime_error("Buffer too small");
            }
            return static_cast<MessageID>(reinterpret_cast<const header_binary*>(_buffer.data())->id);
        }

        const PositionCartesian & get_position_cartesian() const {
            if (get_message_id() != MessageID::MSG_ID_POSITION_CARTESIAN) {
                throw std::runtime_error("Invalid message ID");
            }
            return *reinterpret_cast<const PositionCartesian*>(_buffer.data() + sizeof(header_binary));
        }

        const AttitudeEuler & get_attitude_euler() const {
            if (get_message_id() != MessageID::MSG_ID_ATTITUDE_EULER) {
                throw std::runtime_error("Invalid message ID");
            }
            return *reinterpret_cast<const AttitudeEuler*>(_buffer.data() + sizeof(header_binary));
        }
    
    private:
        enum class ParseState { 
            SYNC, HEADER, BODY
        };

        void handle_one_byte(const char ch, std::function<void()> callback = nullptr) {
            auto header = reinterpret_cast<header_binary*>(_buffer.data());

            _buffer.emplace_back(ch);

            switch (_state)
            {
            case ParseState::SYNC:
                if (_buffer.size() >= 2) {
                    if (header->sync[0] == 0x24 && header->sync[1] == 0x40) {
                        _state = ParseState::HEADER;
                    } else {
                        goto fail_to_parse;
                    }
                }
                break;
            case ParseState::HEADER:/*  */
                if (_buffer.size() >= sizeof(header_binary)) {
                    // header->len is aligned to 4 bytes
                    if (!(header->len & 0x3) && header->len > sizeof(header_binary) && header->len <= MAX_MESSAGE_SIZE) {
                        _state = ParseState::BODY;
                    } else {
                        // invalid message length
                        goto fail_to_parse;
                    }
                }
                break;
            case ParseState::BODY:
                if (_buffer.size() >= header->len) {
                    // check CRC: pass sync and crc16
                    uint16_t crc = crc16((const uint8_t *)&header->id, _buffer.size() - 4);
                    if (crc != header->crc16) {
                        // CRC error
                        goto fail_to_parse;
                    }

                    // message complete
                    if (callback) {
                        callback();
                    }
                    reset();
                }
                break;            
            default: 
            fail_to_parse:
                reset();
            }
        }

        void reset() {
            _buffer.clear();
            _state = ParseState::SYNC;
        }

        ParseState _state {ParseState::SYNC};
        std::vector<char> _buffer;
    };
}

class Component
{
public:
    RCLCPP_SMART_PTR_DEFINITIONS(Component)

    virtual void handle_recv_message(const std::string &bytes) = 0;

    virtual ~Component() = default;
};


class BaseStationComponent : public Component
{
public:
    explicit BaseStationComponent(const rclcpp::Node::SharedPtr &node) : _node(node)
    {
        _bytes_msg.data.reserve(2048);

        _corr_out_pub = _node->create_publisher<std_msgs::msg::ByteMultiArray>("corrections_out", 10);
    }

    void handle_recv_message(const std::string &bytes) override
    {
        auto &data = _bytes_msg.data;
        // check if buffer is full
        if (data.size() + bytes.size() > data.capacity()) {
            RCLCPP_WARN(_node->get_logger(), "Buffer overflow");
            return;
        }

        // append to buffer
        data.insert(data.end(), bytes.begin(), bytes.end());

        // publish when buffer is half full or timeout
        auto pass_time_ms = (_node->now() - _last_pub_time).nanoseconds() / 1e6;
        if (data.size() > 100 || pass_time_ms > 20) {
            _corr_out_pub->publish(_bytes_msg);
            RCLCPP_DEBUG(_node->get_logger(), "Published %zu bytes", data.size());
            data.clear();
            _last_pub_time = _node->now();
        }
    }

private:
    const rclcpp::Node::SharedPtr _node;
    rclcpp::Publisher<std_msgs::msg::ByteMultiArray>::SharedPtr _corr_out_pub;
    std_msgs::msg::ByteMultiArray _bytes_msg {};

    // details
    rclcpp::Time _last_pub_time {_node->now()};
};

#ifndef DEG2RAD
#define DEG2RAD(x) ((x) * 0.01745329251994329575)
#endif

#ifndef RAD2DEG
#define RAD2DEG(x) ((x) * 57.29577951308232087721)
#endif

class RoverComponent : public Component
{
public:
    explicit RoverComponent(const rclcpp::Node::SharedPtr &node, std::function<void(const std::vector<uint8_t> &bytes)> wr_handle) : _node(node), _write_handle(wr_handle)
    {
        const int amc_id = _node->get_parameter("amc_id").as_int();

        _corr_in_sub = _node->create_subscription<std_msgs::msg::ByteMultiArray>("corrections_out", 10,
            [this](std_msgs::msg::ByteMultiArray::SharedPtr msg) {
                RCLCPP_DEBUG(_node->get_logger(), "Received %zu bytes", msg->data.size());
                // send to serial port
                _write_handle(msg->data);
            });
        
        _follower_send_pub = _node->create_publisher<mqsls::msg::FollowerSend>("follower_send" + std::to_string(amc_id), 10);
    }

    void handle_recv_message(const std::string &bytes) override
    {
        _parser.handle_message(bytes, [this]() {
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
                    _follower_send_msg.timestamp = _node->get_clock()->now().nanoseconds() / 1e3; // [us]
                    RCLCPP_INFO(_node->get_logger(), "POSITION_CARTESIAN: %f, %f, %f", _follower_send_msg.position_uav[0], _follower_send_msg.position_uav[1], _follower_send_msg.position_uav[2]);
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

                    RCLCPP_INFO(_node->get_logger(), "ATTITUDE_EULER: heading: %.2f, pitch: %.2f, roll: %.2f", body.Heading, body.Pitch, body.Roll);
                    RCLCPP_INFO(_node->get_logger(), "ATTITUDE_EULER: heading_dot: %.2f, pitch_dot: %.2f, roll_dot: %.2f", body.HeadingDot, body.PitchDot, body.RollDot);

                    _perf_counter.tick();
                    if (_perf_counter.count() % 100 == 0) {
                        _perf_counter.print();
                    }
                    break;
                }
                default:
                    RCLCPP_ERROR(_node->get_logger(), "Unknown message ID: %d", msg_id);
                    break;
            }
        });
    }

private:
    const rclcpp::Node::SharedPtr _node;
    std::function<void(const std::vector<uint8_t> &bytes)> _write_handle;

    rclcpp::Subscription<std_msgs::msg::ByteMultiArray>::SharedPtr _corr_in_sub;
    rclcpp::Publisher<mqsls::msg::FollowerSend>::SharedPtr _follower_send_pub;
    mqsls::msg::FollowerSend _follower_send_msg {};

    // decode
    msg::Parser _parser {};

    PerfCounter _perf_counter {PerfCounterType::INTERVAL, "RTKPlugin"};
};

}