#include <rclcpp/rclcpp.hpp>
#include <gz/math.hh>
#include <gz/msgs.hh>
#include <gz/transport.hh>
#include <gz/msgs/imu.pb.h>
#include <gz/msgs/pose_v.pb.h>
#include <gz/msgs/entity_wrench.pb.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <chrono>
#include <memory>
#include <cfloat>
#include <algorithm>
#include <control_toolbox/pid.hpp>

#include <mqsls/utils.hpp>
#include "data_recorder.hpp"

using namespace std::chrono_literals;
namespace mqsls {

class IdealMqslsControl : public rclcpp::Node {
public:

    /**
     * @brief Construct a new Ideal Mqsls Control object
     */
    IdealMqslsControl(std::chrono::milliseconds control_period = 20ms);

    /**
     * @brief Destroy the Ideal Mqsls Control object
     */
    ~IdealMqslsControl()
    {
        for (auto &sub : _node.SubscribedTopics()) {
            _node.Unsubscribe(sub);
        }
    }

    /**
     * @brief Initialize the node
     * 
     * @return int 
     */
    int init();

    /**
     * @brief rotateQuaternion
     * 
     * @param q_FRD_to_NED 
     * @param q_FLU_to_ENU 
     */
    void rotateQuaternion(gz::math::Quaterniond &q_FRD_to_NED, const gz::math::Quaterniond q_FLU_to_ENU)
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
private:
    // callbacks
    void clockCallback(const gz::msgs::Clock &msg);
    void poseInfoCallback(const gz::msgs::Pose_V &pose);

    // control
    void control_step(const uint64_t dt);

    // data recorder [us]
    void record(const uint64_t timestamp) {
        DataFrame frame;
        frame.timestamp = timestamp;
        frame.load_position = _load_position;
        frame.load_velocity = _load_velocity;
        for (int i = 0; i < 3; i++) {
            frame.uav_position[i] = _uav_position[i];
            frame.uav_velocity[i] = _uav_velocity[i];
        }
        _recorder.push(frame);
    }
    DataRecorder<DataFrame> _recorder{"install/" + utils::realtime_str() + ".csv", 10};

    // publisher for EntityWrench msg
    gz::transport::Node::Publisher _wrench_pub;

    control_toolbox::Pid _pid_z {5, 0.5, 0, 5, -5};

    const std::string _world_name = "mqsls";
    const std::string _payload_name = "ball";
    const std::string _uav_name_prefix = "x500_";
    
    gz::transport::Node _node;
    uint64_t _control_interval; // [us]
    uint64_t _gz_sim_time {0}; // [us]

    // state
    Eigen::Vector3d _load_position;
    Eigen::Vector3d _load_velocity;
    Eigen::Vector3d _uav_position[3];
    Eigen::Vector3d _uav_velocity[3];
};


} // namespace mqsls