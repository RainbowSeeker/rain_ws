#include "fw_formation_control.hpp"

#define M_DEG_TO_RAD 		0.017453292519943295
#define M_RAD_TO_DEG 		57.295779513082323

namespace formation {

FixedwingFormationControl::FixedwingFormationControl(const std::string& node_name, 
                            std::chrono::milliseconds callback_period)
    : Node(node_name)
{
    std::string topic_ns{""};

    if (std::isdigit(node_name.back())) {
        topic_ns = '/' + std::string("px4_")+ node_name.back();
        _uav_id = node_name.back() - '1';

        if (_uav_id < 0 || _uav_id >= 3)
            throw std::invalid_argument("Invalid uav id");
    }

    rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
    auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);
		
    _local_pos_sub = this->create_subscription<px4_msgs::msg::VehicleLocalPosition>(
        topic_ns + "/fmu/out/vehicle_local_position", qos,
        [this](const px4_msgs::msg::VehicleLocalPosition::SharedPtr msg) {
            _local_pos = *msg;
        });
    
    _airspeed_validated_sub = this->create_subscription<px4_msgs::msg::AirspeedValidated>(
        topic_ns + "/fmu/out/airspeed_validated", qos,
        [this](const px4_msgs::msg::AirspeedValidated::SharedPtr msg) {
            _eas2tas = 1.0f; //this is the default value, taken in case of invalid airspeed

            if (std::isfinite(msg->calibrated_airspeed_m_s)
                && std::isfinite(msg->true_airspeed_m_s)
                && (msg->calibrated_airspeed_m_s > 0.0f)) {

                _airspeed = msg->calibrated_airspeed_m_s;

                _eas2tas = constrain(msg->true_airspeed_m_s / msg->calibrated_airspeed_m_s, 0.9f, 2.0f);
            }
        });

    _att_sub = this->create_subscription<px4_msgs::msg::VehicleAttitude>(
        topic_ns + "/fmu/out/vehicle_attitude", qos,
        [this](const px4_msgs::msg::VehicleAttitude::SharedPtr msg) {
            _att = *msg;
        });

    _vehicle_status_sub = this->create_subscription<px4_msgs::msg::VehicleStatus>(
        topic_ns + "/fmu/out/vehicle_status", qos,
        [this](const px4_msgs::msg::VehicleStatus::SharedPtr msg) {
            _vehicle_status = *msg;
        });

    _att_sp_pub = this->create_publisher<px4_msgs::msg::VehicleAttitudeSetpoint>(
        topic_ns + "/fmu/in/vehicle_attitude_setpoint", 10);

    _offboard_control_mode_pub = this->create_publisher<px4_msgs::msg::OffboardControlMode>(
        topic_ns + "/fmu/in/offboard_control_mode", 10);

    _vehicle_command_pub = this->create_publisher<px4_msgs::msg::VehicleCommand>(
        topic_ns + "/fmu/in/vehicle_command", 10);
        
    // Subscribe formation cross
    for (int i = 0; i < 3; i++)
    {
        _form_pos_sub[i] = this->create_subscription<px4_msgs::msg::VehicleLocalPosition>(
        std::string("px4_") + (char)('1' + i) + "/fmu/out/vehicle_local_position", qos,
        [this, i](const px4_msgs::msg::VehicleLocalPosition::SharedPtr msg) {
            _formation_cross.x[i] = msg->x;
            _formation_cross.y[i] = msg->y;
            _formation_cross.h[i] = -msg->z;
            _formation_cross.vn[i] = msg->vx;
            _formation_cross.ve[i] = msg->vy;
            _formation_cross.vd[i] = msg->vz;
            _formation_cross.time_usec[i] = msg->timestamp;
        });
    }
    
    RCLCPP_INFO(this->get_logger(), "Formation node for %s started.", topic_ns.c_str());
    _timer = this->create_wall_timer(callback_period, std::bind(&FixedwingFormationControl::timer_callback, this));
}

FixedwingFormationControl::~FixedwingFormationControl() {
}

void FixedwingFormationControl::timer_callback() {
    // update pilot command & mission data
	pilot_cmd_decode();
	mission_decompose();

	/* only run controller if cross data ready */
	if (formation_preprocess())
	{
		if (_first_ready_time == rclcpp::Time(ROS_ZERO_TIME))
        {
            _first_ready_time = get_clock()->now();
        }

		formation_step();

        running_time = get_clock()->now() - _first_ready_time;
	}
	else
	{
		_first_ready_time = rclcpp::Time(ROS_ZERO_TIME);
	}
}

void FixedwingFormationControl::parameters_update() {
    // FMS parameters
	CONTROL_PARAM.FW_ARSP_MODE = 0;
	CONTROL_PARAM.FW_AIRSPD_MAX = 20.f;
	CONTROL_PARAM.FW_AIRSPD_MIN = 10.f;
	CONTROL_PARAM.FW_AIRSPD_TRIM = 15.f;
	CONTROL_PARAM.FW_AIRSPD_STALL = 7.f;
	CONTROL_PARAM.FW_T_I_GAIN_PIT = 0.1f;
	CONTROL_PARAM.FW_T_I_GAIN_THR = 0.05f;
	CONTROL_PARAM.FW_T_THR_DAMP = 0.1f;
	CONTROL_PARAM.FW_T_SPDWEIGHT = 1.f;
    // CONTROL_PARAM.FW_T_CLMB_R_SP = 3.0f;
    // CONTROL_PARAM.FW_T_SINK_R_SP = 2.0f;
	CONTROL_PARAM.FW_T_CLMB_MAX = 8.f;
	CONTROL_PARAM.FW_T_SINK_MIN = 2.2f;
	CONTROL_PARAM.FW_T_SINK_MAX = 2.7f;
	CONTROL_PARAM.FW_P_LIM_MAX = M_DEG_TO_RAD * 30.f;
	CONTROL_PARAM.FW_P_LIM_MIN = M_DEG_TO_RAD * -15.f;
	CONTROL_PARAM.FW_R_LIM = M_DEG_TO_RAD * 50.f;
	CONTROL_PARAM.FW_T_VERT_ACC = 7.f;
	CONTROL_PARAM.FW_T_RLL2THR = 0.f;
	CONTROL_PARAM.FW_THR_MAX = 0.6f;
	CONTROL_PARAM.FW_THR_MIN = 0.05f;
	CONTROL_PARAM.FW_THR_TRIM = 0.25f;
	CONTROL_PARAM.FW_T_SEB_R_FF = 1.f;
	CONTROL_PARAM.FW_T_PTCH_DAMP = 0.1f;

	// FMS parameters
	FMS_PARAM.FW_AIRSPD_TRIM = CONTROL_PARAM.FW_AIRSPD_TRIM ;
	FMS_PARAM.FW_HEIGHT_TRIM = 100 - FORMATION_PARAM.REL_Z_MATRIX[_uav_id * 3];
	FMS_PARAM.FW_RADIUS_RATIO = 1.0f;
	FMS_PARAM.AIRSPD_P = 1 / 5.0f;
	FMS_PARAM.Z_P = 1 / 5.0f;
	FMS_PARAM.L1 = 50;
	FMS_PARAM.ACCEPT_R = 50;
	FMS_PARAM.LOITER_R = 80;
	FMS_PARAM.ACC_Y_LIM = 5;

	// formation parameters
	FORMATION_PARAM.UAV_ID = _uav_id + 1;
	FORMATION_PARAM.FORM_RADIUS = 1000;
	FORMATION_PARAM.FORM_HGT_KP = 0.5;
	FORMATION_PARAM.FORM_POS_KP = 0.3;
	FORMATION_PARAM.FORM_POUT_LIM = 20;
	FORMATION_PARAM.FORM_VEL_KP = 1.0;
	FORMATION_PARAM.FORM_VOUT_LIM = 2;
	FORMATION_PARAM.FORM_HEAD_KP = 2;
	FORMATION_PARAM.FORM_HOUT_LIM = 0.87266463;
	FORMATION_PARAM.FORM_HRFF_THRE = 20;

	const float hgt_trim = FMS_PARAM.FW_HEIGHT_TRIM;
	// col first			pos1👇	pos2👇	pos3👇
	float form_points[] = { 	0.0f, 	0.0f, 	0.0f,
								200.0f, 200.0f, 200.0f,
								hgt_trim, hgt_trim, hgt_trim};
	for (int i = 0; i < 9; ++i) {
		FORMATION_PARAM.FORM_POINT[i] = form_points[i];
	}

	// col first			pos1👇	pos2👇	pos3👇
	float disband_points[] = {	0.0f, 	00.0f,	00.0f,
								0.0f, 	200.0f, 400.0f,
								hgt_trim, hgt_trim, hgt_trim};
	for (int i = 0; i < 9; ++i) {
		FORMATION_PARAM.DISBAND_POINT[i] = disband_points[i];
	}
}

bool FixedwingFormationControl::formation_preprocess() {
    // check list 1: is flying
    bool is_flying = (_vehicle_status.arming_state == px4_msgs::msg::VehicleStatus::ARMING_STATE_ARMED) &&
            (_vehicle_status.nav_state > px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_MANUAL);

    if (!is_flying)
    {
        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "not flying, nav_state: %d, arming_state: %d", _vehicle_status.nav_state, _vehicle_status.arming_state);
        return false;
    }
#if SIM_MODE == SILSIM
	static rclcpp::Time takeoff_time{get_clock()->now()};
	if (get_clock()->now() - takeoff_time < 8s)
	{
		RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "wait for takeoff.");
		return false;
	}
#endif
	// check list 2: switch to formation mode
	bool is_offboard = _vehicle_status.nav_state == px4_msgs::msg::VehicleStatus::NAVIGATION_STATE_OFFBOARD;
	if (!is_offboard)
	{
		px4_msgs::msg::VehicleCommand cmd{};
		cmd.command = px4_msgs::msg::VehicleCommand::VEHICLE_CMD_DO_SET_MODE;
		cmd.param1 = 1; // base mode
		cmd.param2 = 6;
        cmd.target_system = _vehicle_status.system_id;
		cmd.source_system = _vehicle_status.system_id;
        cmd.target_component = _vehicle_status.component_id;
        cmd.source_component = _vehicle_status.component_id;
        cmd.from_external = true;
        cmd.timestamp = absolute_time();
		_vehicle_command_pub->publish(cmd);

		formation_step();

		RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "Try to switch to formation mode.");
		return false;
	}

	// check list over.
	return true;
}

void FixedwingFormationControl::formation_step() {
    parameters_update();
	/* Run attitude controllers */
	fms_step();
	/* Publish the attitude setpoint for analysis once available */
	publish_attitude_setpoint(_fms_out.att_cmd[0], _fms_out.att_cmd[1], 0, _fms_out.throttle_cmd);
	/* Save data */
	data_save();
}

void FixedwingFormationControl::publish_attitude_setpoint(float roll, float pitch, float yaw, float thrust) {
    px4_msgs::msg::OffboardControlMode ocm{};
	ocm.attitude = true;
	ocm.body_rate = false;
	ocm.timestamp = absolute_time();
	_offboard_control_mode_pub->publish(ocm);

	px4_msgs::msg::VehicleAttitudeSetpoint attitude_setpoint{};
	attitude_setpoint.roll_body = roll;
	attitude_setpoint.pitch_body = pitch;
	attitude_setpoint.yaw_body = yaw;
	attitude_setpoint.thrust_body[0] = thrust;
	attitude_setpoint.timestamp = absolute_time();
	_att_sp_pub->publish(attitude_setpoint);
}

void FixedwingFormationControl::fms_step() {
    // INS_Out
	_fms_in.INS_Out.timestamp = _local_pos.timestamp;
	_fms_in.INS_Out.x_R = _local_pos.x;
	_fms_in.INS_Out.y_R = _local_pos.y;
	_fms_in.INS_Out.h_R = -_local_pos.z;
	_fms_in.INS_Out.vn = _local_pos.vx;
	_fms_in.INS_Out.ve = _local_pos.vy;
	_fms_in.INS_Out.vd = _local_pos.vz;
	_fms_in.INS_Out.airspeed = _airspeed;
    Eigen::Quaterniond q{_att.q[0], _att.q[1], _att.q[2], _att.q[3]};
    double roll, pitch, yaw;
    px4_ros_com::frame_transforms::utils::quaternion::quaternion_to_euler(q, roll, pitch, yaw);
    _fms_in.INS_Out.phi = roll;
    _fms_in.INS_Out.theta = pitch;
    _fms_in.INS_Out.psi = yaw;

	// Formation_Cross
#define ARRAY_LEN(x) (sizeof(x) / sizeof(x[0]))
	for (std::size_t i = 0; i < ARRAY_LEN(_formation_cross.x); i++)
	{
		_fms_in.Formation_Cross.timestamp[i] = _formation_cross.time_usec[i];
		_fms_in.Formation_Cross.x_R[i] = _formation_cross.x[i];
		_fms_in.Formation_Cross.y_R[i] = _formation_cross.y[i];
		_fms_in.Formation_Cross.h_R[i] = _formation_cross.h[i];
		_fms_in.Formation_Cross.vn[i] = _formation_cross.vn[i];
		_fms_in.Formation_Cross.ve[i] = _formation_cross.ve[i];
		_fms_in.Formation_Cross.vd[i] = _formation_cross.vd[i];
	}

	// Run FMS
	_fms_out = _fms.step(&_fms_in);

	auto get_state = [](VehicleState state) -> const char*
	{
		switch (state)
		{
		case VehicleState_None:
			return "None";
		case VehicleState_FormAssemble:
			return "FormAssemble";
		case VehicleState_FormHold:
			return "FormHold";
		case VehicleState_FormMission:
			return "FormMission";
		case VehicleState_FormDisband:
			return "FormDisband";
		case VehicleState_Hold:
			return "Hold";
		default:
			return "UNKNOWN";
		}
	};
    
	auto evaluate_performance = [&](px4_msgs::msg::FormationCross &cross, float rel_x[9], float rel_y[9], float rel_z[9], float psi)
	{
		Eigen::Matrix3f dcm;
		dcm << cosf(psi), sinf(psi), 0.f, -sinf(psi), cosf(psi), 0.f, 0.f, 0.f, 1.f;
		Eigen::Vector3f xyz_err(0, 0, 0);
		for (int i = 0; i < 3; i++)
		{
			Eigen::Vector3f pos_err{cross.x[0] - cross.x[i], cross.y[0] - cross.y[i], cross.h[0] - cross.h[i]};
			pos_err = dcm * pos_err;
			xyz_err(0) += pos_err(0) - rel_x[3 * i];
			xyz_err(1) += pos_err(1) - rel_y[3 * i];
			xyz_err(2) += pos_err(2) - rel_z[3 * i];
		}
		RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000,
			"performance: x_err: %f, y_err: %f, z_err: %f", xyz_err(0), xyz_err(1), xyz_err(2));
	};

	// speed down
	if (_test_phase == PHASE_FORMATION)
	{
		if (_uav_id == 0)
		{
			evaluate_performance(_formation_cross, FORMATION_PARAM.REL_X_MATRIX, FORMATION_PARAM.REL_Y_MATRIX, FORMATION_PARAM.REL_Z_MATRIX, _fms_in.INS_Out.psi);
		}

		RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000,
				"UAV%d state: %s, pose:[%.1f\t%.1f\t%.1f], att:[%.1f\t%.1f\t%.1f]", 
						_uav_id + 1, get_state(_fms_out.FMS_Out.state),
						_fms_in.INS_Out.x_R, _fms_in.INS_Out.y_R, _fms_in.INS_Out.h_R,
						_fms_in.INS_Out.phi * M_RAD_TO_DEG, _fms_in.INS_Out.theta * M_RAD_TO_DEG, _fms_in.INS_Out.psi * M_RAD_TO_DEG);
		
	}
	if (_test_phase == PHASE_SINGLE)
	{
		RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000,
			"setpoint roll: %f, pitch: %f, throttle: %f", M_RAD_TO_DEG * _fms_out.att_cmd[0], M_RAD_TO_DEG * _fms_out.att_cmd[1], _fms_out.throttle_cmd);
	}
}

void FixedwingFormationControl::data_save() {
}

void FixedwingFormationControl::mission_decompose() {
    switch (_fms_out.FMS_Out.state)
	{
	case VehicleState_FormHold:
        if (_fms_in.Mission_Data.timestamp == 0)
        {
            _fms_in.Mission_Data.timestamp = absolute_time();
            _fms_in.Mission_Data.type = 1;
            _fms_in.Mission_Data.valid_items = 2;
            _fms_in.Mission_Data.x[0] = _local_pos.x;
            _fms_in.Mission_Data.y[0] = _local_pos.y;
            _fms_in.Mission_Data.z[0] = _local_pos.z;
            _fms_in.Mission_Data.x[1] = FORMATION_PARAM.FORM_POINT[_uav_id];
            _fms_in.Mission_Data.y[1] = FORMATION_PARAM.FORM_POINT[_uav_id + 3];
            _fms_in.Mission_Data.z[1] = FORMATION_PARAM.FORM_POINT[_uav_id + 6];
        }
		break;
	case VehicleState_FormMission:
		static bool already = false;
		if (_uav_id == 0 && !already)
		{
			already = true;
			_fms_in.Mission_Data.timestamp += 100;
			_fms_in.Mission_Data.type = uint32_t(3);
			_fms_in.Mission_Data.valid_items = 4;
			_fms_in.Mission_Data.x[0] = 0.f;
			_fms_in.Mission_Data.y[0] = 1000.f;
			_fms_in.Mission_Data.z[0] = FMS_PARAM.FW_HEIGHT_TRIM;

			_fms_in.Mission_Data.x[1] = 1000.f;
			_fms_in.Mission_Data.y[1] = 1000.f;
			_fms_in.Mission_Data.z[1] = FMS_PARAM.FW_HEIGHT_TRIM;

			_fms_in.Mission_Data.x[2] = 1000.f;
			_fms_in.Mission_Data.y[2] = 0.f;
			_fms_in.Mission_Data.z[2] = FMS_PARAM.FW_HEIGHT_TRIM;

			_fms_in.Mission_Data.x[3] = 0.f;
			_fms_in.Mission_Data.y[3] = 0.f;
			_fms_in.Mission_Data.z[3] = FMS_PARAM.FW_HEIGHT_TRIM;
		}
		break;
	default:
		break;
	}
}

void FixedwingFormationControl::pilot_cmd_decode() {
    if (_test_phase == PHASE_FORMATION)
	{
		if (running_time < 5min)
		{
			_fms_in.Pilot_Cmd.timestamp = 1;
			_fms_in.Pilot_Cmd.mode = uint32_t(PilotMode_FormAssemble);
		}
		else
		{
			if (_fms_out.FMS_Out.state != VehicleState_FormDisband)
				_fms_in.Pilot_Cmd.timestamp = absolute_time();
			_fms_in.Pilot_Cmd.mode = uint32_t(PilotMode_FormDisband);
		}
	}
}

}; // namespace formation


int main(int argc, const char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <node_name>" << std::endl;
        return 1;
    }

    rclcpp::init(argc, argv);
    
    rclcpp::spin(std::make_shared<formation::FixedwingFormationControl>(argv[1], 40ms));

    rclcpp::shutdown();
    return 0;
}
