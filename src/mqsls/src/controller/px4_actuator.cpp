#include <string>
#include <iostream>
#include <Eigen/Eigen>
#include <rclcpp/rclcpp.hpp>
#include <px4_msgs/msg/vehicle_status.hpp>
#include <px4_msgs/msg/vehicle_command.hpp>
#include <px4_msgs/msg/vehicle_attitude.hpp>
#include <px4_msgs/msg/offboard_control_mode.hpp>
#include <px4_msgs/msg/vehicle_attitude_setpoint.hpp>
#include <mqsls/msg/follower_recv.hpp>

#include "formation/utils.hpp"

namespace mqsls {

#define absolute_time() this->get_clock()->now().nanoseconds() / 1e3

class PX4OutputActuator : public rclcpp::Node
{
public:
    PX4OutputActuator() : Node("px4_output_actuator")
    {
        const std::string topic_ns = "/px4_" + std::to_string(_node_index);

        rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
        auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);

        // Subscriptions
        _follower_recv_sub = this->create_subscription<mqsls::msg::FollowerRecv>(
            "/follower_recv" + std::to_string(_node_index), qos,
            [this](const mqsls::msg::FollowerRecv::SharedPtr msg) {
                apply_force(*msg);
            });


        _vehicle_attitude_sub = this->create_subscription<px4_msgs::msg::VehicleAttitude>(
            topic_ns + "/fmu/out/vehicle_attitude", qos,
            [this](const px4_msgs::msg::VehicleAttitude::SharedPtr msg) {
                _att = *msg;
            });

        _vehicle_status_sub = this->create_subscription<px4_msgs::msg::VehicleStatus>(
            topic_ns + "/fmu/out/vehicle_status", qos,
            [this](const px4_msgs::msg::VehicleStatus::SharedPtr msg) {
                _vehicle_status = *msg;
            });

        // Publishers
        _offboard_control_mode_pub = this->create_publisher<px4_msgs::msg::OffboardControlMode>(
            topic_ns + "/fmu/in/offboard_control_mode", 10);

        _vehicle_attitude_setpoint_pub = this->create_publisher<px4_msgs::msg::VehicleAttitudeSetpoint>(
            topic_ns + "/fmu/in/vehicle_attitude_setpoint", 10);

        _vehicle_command_pub = this->create_publisher<px4_msgs::msg::VehicleCommand>(
            topic_ns + "/fmu/in/vehicle_command", 10);

        RCLCPP_INFO(this->get_logger(), "PX4OutputActuator node for UAV%d has been created.", _node_index);
    }

    void apply_force(const mqsls::msg::FollowerRecv &msg)
    {
        if (!preprocess())
            return;

        if (!_is_init_yaw)
        {
            _init_yaw = utils::quaternion::quaternion_get_yaw(utils::quaternion::array_to_eigen_quat(_att.q));
            _is_init_yaw = true;
        }

        Eigen::Vector3d thrust_sp {msg.force[0], msg.force[1], msg.force[2]};
        Eigen::Quaterniond q_d;
        bodyzToAttitude(-thrust_sp, _init_yaw, q_d);
        
        Eigen::Quaterniond q_att = {_att.q[0], _att.q[1], _att.q[2], _att.q[3]};
        double thrust_project = thrust_sp.dot(q_att.toRotationMatrix() * Eigen::Vector3d(0, 0, -1));
        // convert thrust to normalized thrust. 
        // double thrust_coff = hover_thrust / uav_mass / 9.81;
        publish_attitude_setpoint(q_d, -thrust_project * msg.thrust_coff);
    }
private:
    static void bodyzToAttitude(Eigen::Vector3d body_z, const double yaw_sp, Eigen::Quaterniond &q_d)
    {
        body_z.normalize();

        // vector of desired yaw direction in XY plane, rotated by PI/2
        const Eigen::Vector3d y_C{-sinf(yaw_sp), cosf(yaw_sp), 0.f};

        // desired body_x axis, orthogonal to body_z
        Eigen::Vector3d body_x = y_C.cross(body_z);

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

    // void publish_trajectory_setpoint(Eigen::Vector3d postion)
    // {
    //     px4_msgs::msg::OffboardControlMode ocm{};
    //     ocm.position = true;
    //     ocm.velocity = false;
    //     ocm.acceleration = false;
    //     ocm.attitude = false;
    //     ocm.body_rate = false;
    //     ocm.actuator = false;
    //     ocm.timestamp = absolute_time();
    //     _offboard_control_mode_pub->publish(ocm);

    //     px4_msgs::msg::TrajectorySetpoint setpoint{};
    //     setpoint.position[0] = postion(0);
    //     setpoint.position[1] = postion(1);
    //     setpoint.position[2] = postion(2);
    //     setpoint.yaw = 0;
    //     setpoint.timestamp = absolute_time();
    //     _trajectory_setpoint_pub->publish(setpoint);
    // }

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
        // check list 3: manual control
        if (is_ignore_state(_vehicle_status.nav_state))
        {
            RCLCPP_INFO(this->get_logger(), "UAV%d cannot control, nav_state: %s", _node_index, get_state(_vehicle_status.nav_state));
            return false;
        }

        // check list 4: switch to offboard mode
        bool is_offboard = _vehicle_status.arming_state == px4_msgs::msg::VehicleStatus::ARMING_STATE_ARMED &&
                            _vehicle_status.nav_state == px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_OFFBOARD;
        if (!is_offboard)
        {
            // switch to offboard mode
            publish_vehicle_command(px4_msgs::msg::VehicleCommand::VEHICLE_CMD_DO_SET_MODE, 1, 6);
            publish_vehicle_command(px4_msgs::msg::VehicleCommand::VEHICLE_CMD_COMPONENT_ARM_DISARM, 1);

            RCLCPP_INFO(this->get_logger(), "UAV%d Try to switch to offboard mode.", _node_index);
        }

        return true;
    }

    static bool is_ignore_state(uint8_t state)
    {
        static const uint8_t ignore_state[] = {
            px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_MANUAL,
            px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_LOITER,
            px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_POSCTL,
            px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_LAND,
            px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_RTL,
        };

        for (size_t i = 0; i < sizeof(ignore_state); i++)
        {
            if (state == ignore_state[i])
                return true;
        }
        return false;
    }

    static const char *get_state(uint8_t state)
    {
        static const std::map<uint8_t, std::string> state_str = {
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_MANUAL,         "MANUAL"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_ALTCTL,         "ALTCTL"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_POSCTL,         "POSCTL"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_MISSION,   "AUTO_MISSION"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_LOITER,    "AUTO_LOITER"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_RTL,       "AUTO_RTL"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_ACRO,           "ACRO"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_OFFBOARD,       "OFFBOARD"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_STAB,           "STAB"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_TAKEOFF,   "AUTO_TAKEOFF"},
            {px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_LAND,      "AUTO_LAND"},
        };
        auto it = state_str.find(state);
        if (it != state_str.end())
            return it->second.c_str();
        return "UNKNOWN";
    }

    const int _node_index = this->declare_parameter("amc_id", 1);

    // Subscriptions
    rclcpp::Subscription<mqsls::msg::FollowerRecv>::SharedPtr               _follower_recv_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleAttitude>::SharedPtr         _vehicle_attitude_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleStatus>::SharedPtr           _vehicle_status_sub;

    // Publisher
    rclcpp::Publisher<px4_msgs::msg::OffboardControlMode>::SharedPtr        _offboard_control_mode_pub;
    rclcpp::Publisher<px4_msgs::msg::VehicleAttitudeSetpoint>::SharedPtr    _vehicle_attitude_setpoint_pub;
    rclcpp::Publisher<px4_msgs::msg::VehicleCommand>::SharedPtr             _vehicle_command_pub;

    // store
    px4_msgs::msg::VehicleStatus	    _vehicle_status{};
    px4_msgs::msg::VehicleAttitude	    _att{};
    mqsls::msg::FollowerRecv            _follower_recv{};

    double  _init_yaw = 0;
    bool    _is_init_yaw = false;
};

} // namespace mqsls


int main(int argc, const char** argv) 
{
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<mqsls::PX4OutputActuator>());

    rclcpp::shutdown();

    return 0;
}