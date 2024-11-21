#include <rclcpp/rclcpp.hpp>
#include <gz/math.hh>
#include <gz/msgs.hh>
#include <gz/transport.hh>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <chrono>
#include <memory>
#include <cfloat>
#include <algorithm>
#include <control_toolbox/pid.hpp>
#include <formation/msg/uav_command.hpp>
#include <formation/parameter_manager.hpp>
#include <formation/data_recorder.hpp>
#include <formation/utils.hpp>

namespace mqsls {

struct MqslsDataFrame {
    // timestamp
    uint64_t timestamp;
    // payload state
    Eigen::Vector3d load_position;
    Eigen::Vector3d load_velocity;

    // uav state
    Eigen::Vector3d uav_position[3];
    Eigen::Vector3d uav_velocity[3];

    // operator<<
    friend std::ostream &operator<<(std::ostream &os, const MqslsDataFrame &frame) {
        os << frame.timestamp << ' '
           << frame.load_position.x() << ' ' << frame.load_position.y() << ' ' << frame.load_position.z() << ' '
           << frame.load_velocity.x() << ' ' << frame.load_velocity.y() << ' ' << frame.load_velocity.z() << ' ';
        for (int i = 0; i < 3; i++) {
            os << frame.uav_position[i].x() << ' ' << frame.uav_position[i].y() << ' ' << frame.uav_position[i].z() << ' '
               << frame.uav_velocity[i].x() << ' ' << frame.uav_velocity[i].y() << ' ' << frame.uav_velocity[i].z() << ' ';
        }
        os << '\n';
        return os;
    }

    static std::string header() {
        std::stringstream ss;
        ss  << "timestamp" << ' '
            << "x y z" << ' '
            << "vx vy vz" << ' ';
        for (int i = 1; i <= 3; i++) {
            ss  << "x" << i << " y" << i << " z" << i << ' '
                << "vx" << i << " vy" << i << " vz" << i << ' ';
        }

        return ss.str();
    }
};


class IdealMqslsControl : public rclcpp::Node, public Parameter::ParameterManager
{
public:

    /**
     * @brief Construct a new Ideal Mqsls Control object
     */
    IdealMqslsControl(uint64_t control_period = 20_ms);

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
private:
    
    
    // callbacks
    void clockCallback(const gz::msgs::Clock &msg);
    void poseInfoCallback(const gz::msgs::Pose_V &pose);

    // control
    void control_step(const uint64_t dt);

    // data recorder [us]
    void record(const uint64_t timestamp) {
        MqslsDataFrame frame;
        frame.timestamp = timestamp;
        frame.load_position = _load_position;
        frame.load_velocity = _load_velocity;
        for (int i = 0; i < 3; i++) {
            frame.uav_position[i] = _uav_position[i];
            frame.uav_velocity[i] = _uav_velocity[i];
        }
        _recorder.push(frame);
    }
    DataRecorder<MqslsDataFrame> _recorder{"install/" + utils::nowstr() + ".csv", 10};

    // publisher for EntityWrench msg
    gz::transport::Node::Publisher _wrench_pub;

    // parameters
    Parameter::SharedPtr    _param_world_name {add_parameter("world_name", "mqsls")};

    control_toolbox::Pid _pid_z {5, 0.5, 0, 5, -5};

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