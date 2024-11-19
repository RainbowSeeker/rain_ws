#pragma once

#include <rclcpp/rclcpp.hpp>
#include <cstring>
#include <string>
#include <chrono>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <Eigen/Dense> 
#include <px4_ros_com/frame_transforms.h>
#include <px4_msgs/msg/vehicle_local_position.hpp>
#include <px4_msgs/msg/airspeed_validated.hpp>
#include <px4_msgs/msg/vehicle_attitude.hpp>
#include <px4_msgs/msg/vehicle_status.hpp>
#include <px4_msgs/msg/vehicle_attitude_setpoint.hpp>
#include <px4_msgs/msg/offboard_control_mode.hpp>
#include <px4_msgs/msg/vehicle_command.hpp>
#include <px4_msgs/msg/formation_cross.hpp>

#include <FMS_TECS/FMS.h>

#define SILSIM		1	// Software in the loop simulation
#define REAL		2	// Real flight

#define SIM_MODE	SILSIM

#define PHASE_SINGLE	0
#define PHASE_FORMATION	1
#define PHASE_UNVALID	-1
#define PHASE_DEFAULT	PHASE_FORMATION

#define ROS_ZERO_TIME   0, 0, RCL_ROS_TIME


using namespace std::chrono_literals;

template <typename T>
T constrain(T value, T min_value, T max_value) {
    if (value < min_value) {
        return min_value;
    } else if (value > max_value) {
        return max_value;
    } else {
        return value;
    }
}

namespace formation {

class FixedwingFormationControl : public rclcpp::Node {
public:
    explicit FixedwingFormationControl(const std::string& node_name, 
                            std::chrono::milliseconds callback_period = 20ms);

    ~FixedwingFormationControl();

private:
    inline uint64_t absolute_time() {
        return get_clock()->now().nanoseconds() / 1e3; // [us]
    }

    void timer_callback();
    void parameters_declare();
    void parameters_update();
	bool formation_preprocess();
	void formation_step();
	void publish_attitude_setpoint(float roll, float pitch, float yaw, float thrust);
	void fms_step();
	void data_save();
	void mission_decompose();
	void pilot_cmd_decode();

    rclcpp::TimerBase::SharedPtr _timer;

    // Subscriptions
    rclcpp::Subscription<px4_msgs::msg::VehicleLocalPosition>::SharedPtr    _local_pos_sub;
    rclcpp::Subscription<px4_msgs::msg::AirspeedValidated>::SharedPtr       _airspeed_validated_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleAttitude>::SharedPtr         _att_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleStatus>::SharedPtr           _vehicle_status_sub;
    rclcpp::Subscription<px4_msgs::msg::VehicleLocalPosition>::SharedPtr    _form_pos_sub[3];

    // Publishers
    rclcpp::Publisher<px4_msgs::msg::VehicleAttitudeSetpoint>::SharedPtr    _att_sp_pub;
    rclcpp::Publisher<px4_msgs::msg::OffboardControlMode>::SharedPtr        _offboard_control_mode_pub;
    rclcpp::Publisher<px4_msgs::msg::VehicleCommand>::SharedPtr             _vehicle_command_pub;

    // [m/s] airspeed
	float _airspeed{0.0f};
	float _eas2tas{1.0f};

    px4_msgs::msg::FormationCross		_formation_cross{};
    px4_msgs::msg::VehicleLocalPosition _local_pos{};
	px4_msgs::msg::VehicleAttitude	    _att{};
    px4_msgs::msg::VehicleStatus	    _vehicle_status{};

    int32_t _test_phase{PHASE_DEFAULT};	// determine the running type of uavs. 0: single, 1: formation

    rclcpp::Time _first_ready_time{ROS_ZERO_TIME};
	// rclcpp::Time _last_adhoc_time{ROS_ZERO_TIME};

    int _uav_id{0};
	rclcpp::Duration running_time{0, 0};

    // FMS
	FMS _fms {};
	FMS_In _fms_in {};
	FMS_Out _fms_out {};
};


} // namespace formation

