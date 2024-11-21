#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory.h>
#include <rclcpp/rclcpp.hpp>
#include <string>
#include <chrono>
#include <px4_msgs/msg/vehicle_local_position.hpp>
#include <px4_msgs/msg/vehicle_attitude.hpp>
#include <px4_msgs/msg/vehicle_status.hpp>
#include <px4_msgs/msg/trajectory_setpoint.hpp>
#include <px4_msgs/msg/offboard_control_mode.hpp>
#include <px4_msgs/msg/vehicle_command.hpp>
#include <px4_msgs/msg/formation_cross.hpp>
#include <formation/msg/uav_command.hpp>
#include <formation/msg/uav_status.hpp>
#include <formation/parameter_manager.hpp>
#include <formation/data_recorder.hpp>
#include <formation/utils.hpp>
#include <control_toolbox/pid.hpp>

#include <Eigen/Eigen>
#include "px4_ros_com/frame_transforms.h"



#define SILSIM		1	// Software in the loop simulation
#define REAL		2	// Real flight

#define SIM_MODE	SILSIM

#define PHASE_SINGLE	0
#define PHASE_FORMATION	1

#define ROS_ZERO_TIME   0, 0, RCL_ROS_TIME

using namespace std::chrono_literals;

namespace formation {
class MulticopterFormationControl : public rclcpp::Node, public Parameter::ParameterManager
{

public:
    MulticopterFormationControl(int node_index, std::chrono::milliseconds control_period = 20ms);
    ~MulticopterFormationControl() = default;

private:
    inline uint64_t absolute_time() {
        return get_clock()->now().nanoseconds() / 1e3; // [us]
    }

    inline void update_sta_msg(const std::string &msg) {
        int max_len = std::min(_uav_status.sta_msg.size(), msg.size() + 1);
        std::memcpy(_uav_status.sta_msg.data(), msg.c_str(), max_len);
        _uav_status.sta_msg[max_len - 1] = '\0';
    }

    void timer_callback();
    int  runtime_preprocess();
	int  formation_preprocess();
	void formation_step();
    void formation_enter();
    void formation_exit();
    void fms_step();
    void publish_vehicle_command(uint32_t command, float param1, float param2 = NAN, float param3 = NAN, float param4 = NAN, double param5 = NAN, double param6 = NAN, float param7 = NAN);
	void publish_trajectory_setpoint(float velocity[3], float yaw);
    void handle_command(const formation::msg::UavCommand::SharedPtr msg);

    inline bool uav_is_active() {
        return get_clock()->now() - _last_fmuout_time < 500ms;
    }

    rclcpp::TimerBase::SharedPtr _timer;

    // Subscriptions
    rclcpp::Subscription<px4_msgs::msg::VehicleLocalPosition>::SharedPtr    _local_pos_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleAttitude>::SharedPtr         _att_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleStatus>::SharedPtr           _vehicle_status_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleLocalPosition>::SharedPtr    _form_pos_sub[3];
    // Formation command
    rclcpp::Subscription<formation::msg::UavCommand>::SharedPtr             _command_sub;

    // Publishers
    rclcpp::Publisher<px4_msgs::msg::TrajectorySetpoint>::SharedPtr         _trajectory_setpoint_pub;
    rclcpp::Publisher<px4_msgs::msg::OffboardControlMode>::SharedPtr        _offboard_control_mode_pub;
    rclcpp::Publisher<px4_msgs::msg::VehicleCommand>::SharedPtr             _vehicle_command_pub;
    rclcpp::Publisher<formation::msg::UavStatus>::SharedPtr                 _uav_status_pub;

    // Parameters
    rclcpp::Parameter      _param_test_phase{"test_phase", "formation"};
    rclcpp::Parameter      _param_lasting_time{"lasting_time", 99999}; // [s]
    rclcpp::Parameter      _param_ori_lat{"origin_lat", 47.397742}; // [deg]
    rclcpp::Parameter      _param_ori_lon{"origin_lon", 8.5455940}; // [deg]
    rclcpp::Parameter      _param_ori_alt{"origin_alt", 488.15799}; // [m]
    rclcpp::Parameter      _param_hgt_sp{"mc_hgt_sp", 5.0}; // [m]
    rclcpp::Parameter      _param_hgt_kp{"mc_hgt_kp", 0.5};
    rclcpp::Parameter      _param_hgt_ki{"mc_hgt_ki", 0.05};
    rclcpp::Parameter      _param_yaw_sp{"mc_yaw_sp", 0.0}; // [deg]

    // Deprecated. Parameters are updated in parameter_manager.hpp
    inline void parameters_update()
    {
        _hgt_ctrl.setGains(_param_hgt_kp.as_double(), _param_hgt_ki.as_double(), 0.0, 0.5, -0.5);
    }

    // control_toolbox
    control_toolbox::Pid _hgt_ctrl {_param_hgt_kp.as_double(), _param_hgt_ki.as_double(), 0.0, 0.0, -0.0};

    px4_msgs::msg::FormationCross		_formation_cross{};
    px4_msgs::msg::VehicleLocalPosition _local_pos{};
	px4_msgs::msg::VehicleAttitude	    _att{};
    px4_msgs::msg::VehicleStatus	    _vehicle_status{};
    formation::msg::UavCommand          _command{};
    formation::msg::UavStatus           _uav_status{};
    double _yaw{0.0}; // [rad]
    double _takeoff_hgt{0.0}; // [m]

    int     _test_phase{PHASE_SINGLE};	// determine the running type of uavs. 0: single, 1: formation

    // bool
    bool    _is_stop{true};
    bool    _is_emergency{false};
    bool    _is_landing{false};

    // Time
    rclcpp::Time _first_ready_time{ROS_ZERO_TIME};
    rclcpp::Time _last_fmuout_time{ROS_ZERO_TIME};
    rclcpp::Time _last_cross_time[3]{{ROS_ZERO_TIME}, {ROS_ZERO_TIME}, {ROS_ZERO_TIME}};
    uint64_t _control_interval; // [ns]
    uint64_t _origin_ref_timestamp[3]{UINT64_MAX, UINT64_MAX, UINT64_MAX};
    uint64_t _timer_cycle{0};

    // Duration
	rclcpp::Duration _running_time{0, 0};

    // Formation Related
    int _uav_id{0};
    struct {
        float velocity[3];
        float yaw;
        uint64_t timestamp;
    } _fms_out;
};
} // namespace formation
