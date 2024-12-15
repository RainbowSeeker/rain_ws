#ifndef PX4_COMPONENT_HPP
#define PX4_COMPONENT_HPP

#ifndef COMPONENT_HPP
#error "Include component.hpp instead of px4_component.hpp"
#endif
#include "component.hpp"

#include <px4_msgs/msg/vehicle_status.hpp>
#include <px4_msgs/msg/vehicle_command.hpp>
#include <px4_msgs/msg/vehicle_attitude.hpp>
#include <px4_msgs/msg/trajectory_setpoint.hpp>
#include <px4_msgs/msg/offboard_control_mode.hpp>
#include <px4_msgs/msg/vehicle_local_position.hpp>
#include <px4_msgs/msg/vehicle_attitude_setpoint.hpp>
#include <rclcpp/rclcpp.hpp>
#include <Eigen/Eigen>
#include <iostream>
#include <string>

namespace mqsls {

// Use px4's micro dds to publish 'FollowerSend' msg, 
// so this is a dummy interface class. It does nothing.
class PX4InputSource : public InputSource
{
public:
    PX4InputSource()
    {
        // do nothing
    }

    void register_update_callback(std::function<void(const InputData &)> callback) override
    {
    }
};

#define absolute_time() _node->get_clock()->now().nanoseconds() / 1e3
class PX4OutputActuator : public OutputActuator
{
public:
    PX4OutputActuator(rclcpp::Node *node, int node_index) : _node(node)
    {
        const std::string topic_ns = "/px4_" + std::to_string(node_index);

        rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
        auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);

        // Subscriptions
        _local_pos_sub = _node->create_subscription<px4_msgs::msg::VehicleLocalPosition>(
            topic_ns + "/fmu/out/vehicle_local_position", qos,
            [this](const px4_msgs::msg::VehicleLocalPosition::SharedPtr msg) {
                _local_pos = *msg;
            });

        _vehicle_attitude_sub = _node->create_subscription<px4_msgs::msg::VehicleAttitude>(
            topic_ns + "/fmu/out/vehicle_attitude", qos,
            [this](const px4_msgs::msg::VehicleAttitude::SharedPtr msg) {
                _att = *msg;
            });

        _vehicle_status_sub = _node->create_subscription<px4_msgs::msg::VehicleStatus>(
            topic_ns + "/fmu/out/vehicle_status", qos,
            [this](const px4_msgs::msg::VehicleStatus::SharedPtr msg) {
                _vehicle_status = *msg;
            });

        // Publishers
        _offboard_control_mode_pub = _node->create_publisher<px4_msgs::msg::OffboardControlMode>(
            topic_ns + "/fmu/in/offboard_control_mode", 10);

        _vehicle_attitude_setpoint_pub = _node->create_publisher<px4_msgs::msg::VehicleAttitudeSetpoint>(
            topic_ns + "/fmu/in/vehicle_attitude_setpoint", 10);

        _trajectory_setpoint_pub = _node->create_publisher<px4_msgs::msg::TrajectorySetpoint>(
            topic_ns + "/fmu/in/trajectory_setpoint", 10);

        _vehicle_command_pub = _node->create_publisher<px4_msgs::msg::VehicleCommand>(
            topic_ns + "/fmu/in/vehicle_command", 10);

        // Parameter declare
        if (!_node->has_parameter("hover_thrust")) _node->declare_parameter("hover_thrust", 0.5);
        if (!_node->has_parameter("uav_mass")) _node->declare_parameter("uav_mass", 1.0);
    }

    void apply(const OutputData &data) override
    {
        if (!preprocess())
            return;

        // Get parameters
        const double hover_thrust = _node->get_parameter("hover_thrust").as_double();
        const double uav_mass = _node->get_parameter("uav_mass").as_double();

        Eigen::Vector3d thrust_sp {data.msg.force[0], data.msg.force[1], data.msg.force[2]};
        Eigen::Quaterniond q_d;
        bodyzToAttitude(-thrust_sp, q_d);
        
        Eigen::Quaterniond q_att = {_att.q[0], _att.q[1], _att.q[2], _att.q[3]};
        double thrust_project = thrust_sp.dot(q_att.toRotationMatrix() * Eigen::Vector3d(0, 0, -1));
        // convert thrust to normalized thrust. 
        double thrust_coff = hover_thrust / uav_mass / 9.81;
        publish_attitude_setpoint(q_d, -thrust_project * thrust_coff);
    }
private:
    static void bodyzToAttitude(Eigen::Vector3d body_z, Eigen::Quaterniond &q_d)
    {
        body_z.normalize();

        const Eigen::Vector3d y{0, 1, 0};

        // desired body_x axis, orthogonal to body_z
        Eigen::Vector3d body_x = y.cross(body_z);

        // keep nose to front while inverted upside down
        if (body_z(2) < 0.0f) {
            body_x = -body_x;
        }

        body_x.normalize();

        // desired body_y axis
        const Eigen::Vector3d body_y = body_z.cross(body_x);

        Eigen::Matrix3d R_sp;

        // fill rotation matrix
        for (int i = 0; i < 3; i++) {
            R_sp(i, 0) = body_x(i);
            R_sp(i, 1) = body_y(i);
            R_sp(i, 2) = body_z(i);
        }

        // copy quaternion setpoint to attitude setpoint topic
        q_d = Eigen::Quaterniond(R_sp);
    }

    void publish_attitude_setpoint(Eigen::Quaterniond q_d, double thrust)
    {
        px4_msgs::msg::OffboardControlMode ocm{};
        ocm.position = false;
        ocm.velocity = false;
        ocm.acceleration = false;
        ocm.attitude = true;
        ocm.body_rate = false;
        ocm.actuator = false;
        ocm.timestamp = absolute_time();
        _offboard_control_mode_pub->publish(ocm);

        px4_msgs::msg::VehicleAttitudeSetpoint att_sp{};
        att_sp.q_d[0] = q_d.w();
        att_sp.q_d[1] = q_d.x();
        att_sp.q_d[2] = q_d.y();
        att_sp.q_d[3] = q_d.z();
        att_sp.thrust_body[0] = 0.0f;
        att_sp.thrust_body[1] = 0.0f;
        att_sp.thrust_body[2] = thrust; // only used in multicopter
        att_sp.yaw_sp_move_rate = NAN;
        att_sp.reset_integral = false;
        att_sp.fw_control_yaw_wheel = false;
        att_sp.timestamp = absolute_time();
        _vehicle_attitude_setpoint_pub->publish(att_sp);
    }

    void publish_trajectory_setpoint(Eigen::Vector3d postion)
    {
        px4_msgs::msg::OffboardControlMode ocm{};
        ocm.position = true;
        ocm.velocity = false;
        ocm.acceleration = false;
        ocm.attitude = false;
        ocm.body_rate = false;
        ocm.actuator = false;
        ocm.timestamp = absolute_time();
        _offboard_control_mode_pub->publish(ocm);

        px4_msgs::msg::TrajectorySetpoint setpoint{};
        setpoint.position[0] = postion(0);
        setpoint.position[1] = postion(1);
        setpoint.position[2] = postion(2);
        setpoint.yaw = 0;
        setpoint.timestamp = absolute_time();
        _trajectory_setpoint_pub->publish(setpoint);
    }

    void publish_vehicle_command(uint32_t command, float param1, float param2 = NAN, float param3 = NAN, float param4 = NAN, double param5 = NAN, double param6 = NAN, float param7 = NAN)
    {
        px4_msgs::msg::VehicleCommand cmd{};
        cmd.command = command;
        cmd.param1 = param1;
        cmd.param2 = param2;
        cmd.param3 = param3;
        cmd.param4 = param4;
        cmd.param5 = param5;
        cmd.param6 = param6;
        cmd.param7 = param7;
        cmd.target_system = _vehicle_status.system_id;
        cmd.source_system = _vehicle_status.system_id;
        cmd.target_component = _vehicle_status.component_id;
        cmd.source_component = _vehicle_status.component_id;
        cmd.from_external = true;
        cmd.timestamp = absolute_time();
        _vehicle_command_pub->publish(cmd);
    }

    bool preprocess()
    {
        // check list 4: switch to offboard mode
        bool is_offboard = _vehicle_status.arming_state == px4_msgs::msg::VehicleStatus::ARMING_STATE_ARMED &&
                            _vehicle_status.nav_state == px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_OFFBOARD;
        if (!is_offboard)
        {
            // switch to offboard mode
            publish_vehicle_command(px4_msgs::msg::VehicleCommand::VEHICLE_CMD_DO_SET_MODE, 1, 6);
            publish_vehicle_command(px4_msgs::msg::VehicleCommand::VEHICLE_CMD_COMPONENT_ARM_DISARM, 1);

            // takeoff
            publish_attitude_setpoint(Eigen::Quaterniond::Identity(), -_node->get_parameter("hover_thrust").as_double());

            RCLCPP_INFO_THROTTLE(_node->get_logger(), *_node->get_clock(), 1000, "Try to switch to offboard mode.");
            return false;
        }

        return true;
    }

    rclcpp::Node *_node;

    // Subscriptions
    rclcpp::Subscription<px4_msgs::msg::VehicleLocalPosition>::SharedPtr    _local_pos_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleAttitude>::SharedPtr         _vehicle_attitude_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleStatus>::SharedPtr           _vehicle_status_sub;

    // Publisher
    rclcpp::Publisher<px4_msgs::msg::OffboardControlMode>::SharedPtr        _offboard_control_mode_pub;
    rclcpp::Publisher<px4_msgs::msg::VehicleAttitudeSetpoint>::SharedPtr    _vehicle_attitude_setpoint_pub;
    rclcpp::Publisher<px4_msgs::msg::TrajectorySetpoint>::SharedPtr         _trajectory_setpoint_pub;
    rclcpp::Publisher<px4_msgs::msg::VehicleCommand>::SharedPtr             _vehicle_command_pub;

    // store
    px4_msgs::msg::VehicleStatus	    _vehicle_status{};
    px4_msgs::msg::VehicleLocalPosition _local_pos{};
    px4_msgs::msg::VehicleAttitude	    _att{};
};


class PX4EventHandler : public EventHandler
{
public:
    PX4EventHandler(rclcpp::Node *node) : _node(node)
    {
    }

    void register_periodic_callback(std::function<void(uint64_t)> callback, uint64_t period_us) override
    {
        if (_timer)
        {
            _timer->cancel();
            RCLCPP_ERROR(_node->get_logger(), "Timer already exists. Cancel it.");
        }
        _callback = callback;
        _timer = _node->create_wall_timer(std::chrono::microseconds(period_us), std::bind(&PX4EventHandler::timer_callback, this));
    }
private:
    void timer_callback()
    {
        _callback(absolute_time());
    }

    rclcpp::Node *_node;

    // timer
    rclcpp::TimerBase::SharedPtr _timer;
    std::function<void(uint64_t)> _callback;
};


} // namespace mqsls

#endif // !PX4_COMPONENT_HPP