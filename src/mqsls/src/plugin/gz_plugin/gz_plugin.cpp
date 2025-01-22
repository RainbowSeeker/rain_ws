#include <string>
#include <iostream>
#include <gz/msgs.hh>
#include <gz/math.hh>
#include <gz/transport.hh>
#include <rclcpp/rclcpp.hpp>
#include <Eigen/Eigen>
#include <mqsls/msg/follower_send.hpp>

namespace mqsls {

class GzSensor : public rclcpp::Node
{
public:
    GzSensor() : rclcpp::Node("GzSensor")
    {
        const std::string world_name = this->declare_parameter("world_name", "default");

        // pose: /world/$WORLD/pose/info
        std::string pose_topic = "/world/" + world_name + "/pose/info";
        if (!_node.Subscribe(pose_topic, &GzSensor::poseInfoCallback, this)) {
            throw std::runtime_error("Failed to subscribe to " + pose_topic);
        }

        // Publisher
        for (int i = 0; i < 3; i++) {
            _follower_send_pub[i] = this->create_publisher<mqsls::msg::FollowerSend>(
                "follower_send" + std::to_string(i + 1), 10);
        }

        _timer = this->create_wall_timer(std::chrono::milliseconds(5), std::bind(&GzSensor::run, this));
    }

private:
    struct pos_state
    {
        uint64_t last_time;
        Eigen::Vector3d position;
        Eigen::Vector3d velocity;
    };

    /**
     * @brief rotateQuaternion
     * 
     * @param q_FRD_to_NED 
     * @param q_FLU_to_ENU 
     */
    static void rotateQuaternion(gz::math::Quaterniond &q_FRD_to_NED, const gz::math::Quaterniond q_FLU_to_ENU)
    {
        // FLU (ROS) to FRD (PX4) static rotation
        static const auto q_FLU_to_FRD = gz::math::Quaterniond(0, 1, 0, 0);

        /**
         * @brief Quaternion for rotation between ENU and NED frames
         *
         * NED to ENU: +PI/2 rotation about Z (Down) followed by a +PI rotation around X (old North/new East)
         * ENU to NED: +PI/2 rotation about Z (Up) followed by a +PI rotation about X (old East/new North)
         * This rotation is symmetric, so q_ENU_to_NED == q_NED_to_ENU.
         */
        static const auto q_ENU_to_NED = gz::math::Quaterniond(0, 0.70711, 0.70711, 0);

        // final rotation composition
        q_FRD_to_NED = q_ENU_to_NED * q_FLU_to_ENU * q_FLU_to_FRD.Inverse();
    }

    void poseInfoCallback(const gz::msgs::Pose_V &pose)
    {
    #define get_ned_motion(ref_last_time, ref_postion, ref_velocity) \
    { \
        const double dt = std::clamp((current_time - ref_last_time) * 1e-6, 0.001, 0.1); \
        ref_last_time = current_time; \
        gz::msgs::Vector3d pose_position = pose.pose(i).position(); \
        gz::msgs::Quaternion pose_orientation = pose.pose(i).orientation(); \
        gz::math::Quaterniond q_gr = gz::math::Quaterniond( \
                                pose_orientation.w(), \
                                pose_orientation.x(), \
                                pose_orientation.y(), \
                                pose_orientation.z()); \
        gz::math::Quaterniond q_nb; \
        rotateQuaternion(q_nb, q_gr); \
        const Eigen::Vector3d position{pose_position.y(), pose_position.x(), -pose_position.z()}; \
        const Eigen::Vector3d velocity{(position - ref_postion) / dt}; \
        const Eigen::Vector3d acceleration{(velocity - ref_velocity) / dt}; \
        const Eigen::Quaterniond att(q_nb.W(), q_nb.X(), q_nb.Y(), q_nb.Z()); \
        ref_postion = position; \
        ref_velocity = velocity; \
        }
        //const Eigen::Vector3d euler = att.toRotationMatrix().eulerAngles(2, 1, 0);
        //const Eigen::Vector3d angular_rate = (euler - ref_euler) / dt;
        //ref_euler = euler;
        //ref_angular_rate = angular_rate;
        static const Eigen::Vector3d max_velocity{10, 10, 10}; 
        static const Eigen::Vector3d min_velocity{-10, -10, -10};
        const uint64_t current_time = pose.header().stamp().sec() * 1e6 + pose.header().stamp().nsec() / 1e3;

        for (int i = 0; i < pose.pose_size(); i++) {
            if (pose.pose(i).name() == _payload_name) {
                
                get_ned_motion(_load_state.last_time, _load_state.position, _load_state.velocity);

                // saturation
                _load_state.velocity = _load_state.velocity.cwiseMax(min_velocity).cwiseMin(max_velocity);

            }
            else if (pose.pose(i).name().substr(0, _uav_name_prefix.size()) == _uav_name_prefix) {
                const int uav_index = std::stoi(pose.pose(i).name().substr(_uav_name_prefix.size())) - 1;
                // idx \in [0, 1, 2]
                if (uav_index < 0 || uav_index >= 3) {
                    std::cout << "Invalid UAV index: " << uav_index << std::endl;
                }

                get_ned_motion(_uav_state[uav_index].last_time, _uav_state[uav_index].position, _uav_state[uav_index].velocity);

                // saturation
                _uav_state[uav_index].velocity = _uav_state[uav_index].velocity.cwiseMax(min_velocity).cwiseMin(max_velocity);
            }
        }

    #undef get_ned_motion
    }

    void run()
    {
        for (int i = 0; i < 3; i++) {
            mqsls::msg::FollowerSend msg;
            msg.timestamp = _load_state.last_time;
            for (int j = 0; j < 3; j++) {
                msg.position_load[j] = _load_state.position[j];
                msg.velocity_load[j] = _load_state.velocity[j];
                msg.position_uav[j] = _uav_state[i].position[j];
                msg.velocity_uav[j] = _uav_state[i].velocity[j];
            }
            _follower_send_pub[i]->publish(msg);
        }
    }

    // Gz
    gz::transport::Node _node;
    const std::string _payload_name = "payload";
    const std::string _uav_name_prefix = "x500_";

    // ros2
    rclcpp::TimerBase::SharedPtr _timer;
    rclcpp::Publisher<mqsls::msg::FollowerSend>::SharedPtr _follower_send_pub[3];

    // detail
    pos_state _load_state, _uav_state[3];
};

} // namespace mqsls

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<mqsls::GzSensor>());
    rclcpp::shutdown();
    return 0;
}