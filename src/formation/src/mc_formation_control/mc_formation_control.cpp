#include "mc_formation_control.hpp"
#include "formation/utils.hpp"
#include "px4_ros_com/frame_transforms.h"
#include <variant>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

using namespace Eigen;
using namespace px4_msgs::msg;

namespace formation {

// static
static const int RESULT_SUCCESS     = 0;
static const int RESULT_WAITING     = 1;
static const int RESULT_FAILED      = -1;
static const uint8_t ignore_state[] = {
    px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_MANUAL,
    px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_LOITER,
    px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_POSCTL,
    px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_LAND,
    px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_AUTO_RTL,
};
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

static const char *get_state(uint8_t state)
{
    auto it = state_str.find(state);
    if (it != state_str.end())
        return it->second.c_str();
    return "UNKNOWN";
}

static bool is_ignore_state(uint8_t state)
{
    for (size_t i = 0; i < sizeof(ignore_state); i++)
    {
        if (state == ignore_state[i])
            return true;
    }
    return false;
}

MulticopterFormationControl::MulticopterFormationControl(int node_index, std::chrono::milliseconds control_period) :
    Node("amc_" + std::to_string(node_index)), Parameter::ParameterManager(this) ,_control_interval(control_period / 1ms * 1e6) // [ns]
{
    std::string topic_ns{""};

    if (node_index)
    {
        _uav_id = node_index - 1;
        if (_uav_id < 0 || _uav_id >= 3)
            throw std::invalid_argument("Invalid node index");

        topic_ns = "/px4_" + std::to_string(node_index);
    }

    rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
    auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);

    _local_pos_sub = this->create_subscription<VehicleLocalPosition>(
        topic_ns + "/fmu/out/vehicle_local_position", qos,
        [this](const VehicleLocalPosition::SharedPtr msg) {
            _local_pos = *msg;
            _last_fmuout_time = get_clock()->now();
        });

    _att_sub = this->create_subscription<VehicleAttitude>(
        topic_ns + "/fmu/out/vehicle_attitude", qos,
        [this](const VehicleAttitude::SharedPtr msg) {
            using namespace px4_ros_com::frame_transforms::utils::quaternion;
            _att = *msg;
            _yaw = quaternion_get_yaw(array_to_eigen_quat(_att.q));
        });

    _vehicle_status_sub = this->create_subscription<VehicleStatus>(
        topic_ns + "/fmu/out/vehicle_status", qos,
        [this](const VehicleStatus::SharedPtr msg) {
            _vehicle_status = *msg;
        });

    _trajectory_setpoint_pub = this->create_publisher<TrajectorySetpoint>(
        topic_ns + "/fmu/in/trajectory_setpoint", 10);

    _offboard_control_mode_pub = this->create_publisher<OffboardControlMode>(
        topic_ns + "/fmu/in/offboard_control_mode", 10);

    _vehicle_command_pub = this->create_publisher<VehicleCommand>(
        topic_ns + "/fmu/in/vehicle_command", 10);

    // Subscribe formation cross
    _test_phase = _param_test_phase->as_string() == "formation" ? PHASE_FORMATION : PHASE_SINGLE;

    if (_test_phase == PHASE_FORMATION)
    {
        for (int i = 0; i < 3; i++)
        {
            _form_pos_sub[i] = this->create_subscription<VehicleLocalPosition>(
            std::string("px4_") + (char)('1' + i) + "/fmu/out/vehicle_local_position", qos,
            [this, i](const VehicleLocalPosition::SharedPtr msg) {
                _formation_cross.x[i] = msg->x;
                _formation_cross.y[i] = msg->y;
                _formation_cross.h[i] = -msg->z;
                _formation_cross.vn[i] = msg->vx;
                _formation_cross.ve[i] = msg->vy;
                _formation_cross.vd[i] = msg->vz;
                _formation_cross.time_usec[i] = msg->timestamp;
                _last_cross_time[i] = get_clock()->now();

                // check ekf origin is same
                if (_origin_ref_timestamp[i] != msg->timestamp)
                {
                    _origin_ref_timestamp[i] = msg->timestamp;
                    _is_same_origin[i] = false;
                }

                if (!_is_same_origin[i] && std::isfinite(msg->ref_lat) && std::abs(msg->ref_lat - _param_ori_lat->as_double()) < 1e-7 &&
                    std::isfinite(msg->ref_lon) && std::abs(msg->ref_lon - _param_ori_lon->as_double()) < 1e-7 &&
                    std::isfinite(msg->ref_alt) && std::abs(msg->ref_alt - (float)_param_ori_alt->as_double()) < 1e-3)
                {
                    _is_same_origin[i] = true;
                }
            });
        }
    }

    _command_sub = this->create_subscription<formation::msg::UavCommand>(
        topic_ns + "/fmu/in/uav_command", 10,
        std::bind(&MulticopterFormationControl::handle_command, this, std::placeholders::_1));

    _uav_status_pub = this->create_publisher<formation::msg::UavStatus>(
        topic_ns + "/fmu/out/uav_status", 10);

    RCLCPP_INFO(this->get_logger(), "Formation node for %s started. Phase: %s", topic_ns.c_str(), _param_test_phase->as_string().c_str());
    _timer = this->create_wall_timer(control_period, std::bind(&MulticopterFormationControl::timer_callback, this));
}


void MulticopterFormationControl::timer_callback()
{
    if (runtime_preprocess())
        return;

	/* only run controller if cross data ready */
    switch (formation_preprocess())
    {
    case RESULT_SUCCESS:
        formation_enter();

        formation_step();
        break;
    case RESULT_WAITING:
        // do nothing
        break;
    case RESULT_FAILED:
    default:
        formation_exit();
        break;
    }
}

int MulticopterFormationControl::runtime_preprocess()
{
    _timer_cycle += 1;
    // publish uav status
    if (_timer_cycle % 2 == 0)
    {
        _uav_status.timestamp = absolute_time();
        _uav_status.stage = uav_is_active() && _vehicle_status.nav_state == VehicleStatus::NAVIGATION_STATE_OFFBOARD ? 1 : 0;
        _uav_status_pub->publish(_uav_status);
    }

    // is over time ?
    if (_running_time > 1s * _param_lasting_time->as_int())
    {
        if (_vehicle_status.nav_state == VehicleStatus::NAVIGATION_STATE_AUTO_LAND)
        {
            _is_landing = true;
        }

        if (!_is_landing)
        {
            // switch to land mode
            publish_vehicle_command(VehicleCommand::VEHICLE_CMD_DO_SET_MODE, 1, 4, 6);
            RCLCPP_WARN_THROTTLE(this->get_logger(), *get_clock(), 1000, "Timeout: Try to switch to land mode.");
        }
        update_sta_msg("Landing.");
        return RESULT_WAITING;
    }

    return RESULT_SUCCESS;
}

int MulticopterFormationControl::formation_preprocess()
{
#if SIM_MODE == REAL
    // check list 0: is ready
    if (_is_stop)
    {
        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "Wait for GCS.");
        update_sta_msg("Wait for GCS.");
        return RESULT_FAILED;
    }
#endif

    // check list 1: px4 data ready
    if (!uav_is_active())
    {
        RCLCPP_WARN_THROTTLE(this->get_logger(), *get_clock(), 1000, "Fmu data lost.");
        update_sta_msg("Fmu data lost.");
        return RESULT_FAILED;
    }

    // check list 2: cross data ready
    if (_test_phase == PHASE_FORMATION)
    {
        for (int i = 0; i < 3; i++)
        {
            // check origin is same
            if (!_is_same_origin[i])
            {
                RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "Ekf origin is not same with [%d].", i + 1);
                update_sta_msg("Ekf ori isn't same.");
                return RESULT_FAILED;
            }

            if (get_clock()->now() - _last_cross_time[i] > 1s)
            {
                RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "Wait for cross data from id: %d.", i + 1);
                update_sta_msg("Wait for cross data.");
                return RESULT_FAILED;
            }
        }
    }

    // check list 3: manual control
    if (is_ignore_state(_vehicle_status.nav_state))
    {
        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "cannot control, nav_state: %s", get_state(_vehicle_status.nav_state));
        update_sta_msg("Cannot control.");
        return RESULT_WAITING;
    }

	// check list 4: switch to formation mode
	bool is_offboard = _vehicle_status.arming_state == VehicleStatus::ARMING_STATE_ARMED &&
                        _vehicle_status.nav_state == VehicleStatus::NAVIGATION_STATE_OFFBOARD;
	if (!is_offboard)
	{
		// switch to offboard mode
        publish_vehicle_command(VehicleCommand::VEHICLE_CMD_DO_SET_MODE, 1, 6);
        publish_vehicle_command(VehicleCommand::VEHICLE_CMD_COMPONENT_ARM_DISARM, 1);

        float velocity[3] = {0.0f, 0.0f, -0.5f};
        publish_trajectory_setpoint(velocity, _yaw);

		RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "Try to switch to formation mode.");
		update_sta_msg("switch to formation.");
		return RESULT_WAITING;
	}

	// check list over.
	update_sta_msg("Normal.");
	return RESULT_SUCCESS;
}

void MulticopterFormationControl::formation_enter()
{
    if (_first_ready_time == rclcpp::Time(ROS_ZERO_TIME))
    {
        _first_ready_time = get_clock()->now();
        _takeoff_hgt = -_local_pos.z;
    }

    _running_time = get_clock()->now() - _first_ready_time;
}

void MulticopterFormationControl::formation_exit()
{
    // when uav firstly get ready and exit, enter state of emergency.
    if (_first_ready_time != rclcpp::Time(ROS_ZERO_TIME))
    {
        _is_emergency = true;
        _first_ready_time = rclcpp::Time(ROS_ZERO_TIME);
    }

    if (_is_emergency && uav_is_active())
    {
        if (_vehicle_status.nav_state == VehicleStatus::NAVIGATION_STATE_OFFBOARD
            || !is_ignore_state(_vehicle_status.nav_state))
        {
            // switch to hold mode
            publish_vehicle_command(VehicleCommand::VEHICLE_CMD_DO_SET_MODE, 1, 4, 3);
            RCLCPP_WARN_THROTTLE(this->get_logger(), *get_clock(), 1000,
                    "emergency: Try to switch to hold mode.");
        } else
        {
            _is_emergency = false;
            RCLCPP_INFO(this->get_logger(), "emergency: Succeed in escaping from danger.");
        }
    }
}

void MulticopterFormationControl::formation_step()
{
    parameters_update();
	/* Run attitude controllers */
	fms_step();
	/* Publish the attitude setpoint for analysis once available */
	publish_trajectory_setpoint(_fms_out.velocity, _fms_out.yaw);
}

void MulticopterFormationControl::fms_step()
{
    // control interval [ns]
    const uint64_t dt = _control_interval;

    // height hold
    double vh_cmd = _hgt_ctrl.computeCommand(_takeoff_hgt + _param_hgt_sp->as_double() + _local_pos.z, dt);

    // yaw hold
    const double yaw_cmd = math::radians(_param_yaw_sp->as_double());

    // formation control
    Vector3d pos_err{0, 0, 0};
    Vector3d vel_sp {0.25, 0, 0};

    if (_test_phase == PHASE_SINGLE)
    {
        // follow yaw.
        Matrix3d Rz;
        Rz << cos(_yaw), -sin(_yaw),  0,
              sin(_yaw),  cos(_yaw),  0,
              0,                  0,  1;

        vel_sp = Rz * vel_sp;
        vel_sp[2] += -vh_cmd;
    }
    else if (_test_phase == PHASE_FORMATION)
    {
        Matrix3d rel_x;
        rel_x << 0.0, -2.0, -2.0,
                2.0, 0.0, 0.0,
                2.0, 0.0, 0.0;
        Matrix3d rel_y;
        rel_y << 0.0, -2.0, 2.0,
                2.0, 0.0, 4.0,
                -2.0, -4.0, 0.0;
        Matrix3d rel_z = Matrix3d::Zero();
        for (int i = 0; i < 3; i++)
        {
            pos_err += Vector3d{_formation_cross.x[i] + rel_x(i, _uav_id), _formation_cross.y[i] + rel_y(i, _uav_id), -_formation_cross.h[i] + rel_z(i, _uav_id)} / 3.0;
            // vel_sp  += Vector3d{_formation_cross.vn[i], _formation_cross.ve[i], _formation_cross.vd[i]} / 3.0;
        }
        pos_err -= Vector3d{_local_pos.x, _local_pos.y, _local_pos.z};
        vel_sp += pos_err * 0.5;
        vel_sp[2] += -vh_cmd;
    }

    // output limit
    vel_sp[0] = math::constrain(vel_sp[0], -0.5, 0.5);
    vel_sp[1] = math::constrain(vel_sp[1], -0.5, 0.5);
    vel_sp[2] = math::constrain(vel_sp[2], -0.5, 0.5);

    _fms_out.velocity[0] = (float)vel_sp[0];
    _fms_out.velocity[1] = (float)vel_sp[1];
    _fms_out.velocity[2] = (float)vel_sp[2];
    _fms_out.yaw = (float)yaw_cmd;

    if (_test_phase == PHASE_SINGLE)
    {
        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000,
            "command vn: %f, ve: %f, vd: %f, yaw_cmd: %f deg", vel_sp[0], vel_sp[1], vel_sp[2], math::degrees(yaw_cmd));
    }
    else if (_test_phase == PHASE_FORMATION)
    {
        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000,
            "formation: x_err: %f, y_err: %f, z_err: %f", pos_err(0), pos_err(1), pos_err(2));
    }
}

void MulticopterFormationControl::publish_vehicle_command(uint32_t command, float param1, float param2, float param3, float param4, double param5, double param6, float param7)
{
    VehicleCommand cmd{};
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

void MulticopterFormationControl::publish_trajectory_setpoint(float velocity[3], float yaw)
{
    OffboardControlMode ocm{};
	ocm.position = false;
	ocm.velocity = true;
	ocm.acceleration = false;
	ocm.attitude = false;
	ocm.body_rate = false;
	ocm.actuator = false;
	ocm.timestamp = absolute_time();
	_offboard_control_mode_pub->publish(ocm);

	TrajectorySetpoint setpoint{};
    setpoint.position = {NAN, NAN, NAN};
    setpoint.velocity = {velocity[0], velocity[1], velocity[2]};
    setpoint.yaw = yaw;
    setpoint.timestamp = absolute_time();
	_trajectory_setpoint_pub->publish(setpoint);
}

void MulticopterFormationControl::handle_command(const formation::msg::UavCommand::SharedPtr msg)
{
    using namespace formation::msg;

    switch (msg->command)
    {
    case UavCommand::UAV_CMD_FORM_START:
        _is_stop = false;
        break;
    case UavCommand::UAV_CMD_FORM_STOP:
        _is_stop = true;
        break;
    default:
        break;
    }

    RCLCPP_INFO(this->get_logger(), "Received command: %d, param: %.3f %.3f", msg->command, msg->param1, msg->param2);
}

} // namespace formation

int main(int argc, const char** argv)
{
    int node_index = argc > 1 ? std::atoi(argv[1]) : 0;

    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<formation::MulticopterFormationControl>(node_index, 40ms));

    rclcpp::shutdown();
    return 0;
}
