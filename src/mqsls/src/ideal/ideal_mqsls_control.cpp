#include "ideal_mqsls_control.hpp"

namespace mqsls {

IdealMqslsControl::IdealMqslsControl(uint64_t control_period) : 
    rclcpp::Node("ideal_mqsls_control"), Parameter::ParameterManager(this), _control_interval(control_period)
{
    init();
}

int IdealMqslsControl::init()
{
    const std::string world_name = _param_world_name->as_string();
    
    // clock
    std::string clock_topic = "/world/" + world_name + "/clock";
    if (!_node.Subscribe(clock_topic, &IdealMqslsControl::clockCallback, this)) {
        RCLCPP_ERROR(this->get_logger(), "Failed to subscribe to %s", clock_topic.c_str());
        return -1;
    }

    // pose: /world/$WORLD/pose/info
    std::string pose_topic = "/world/" + world_name + "/pose/info";
    if (!_node.Subscribe(pose_topic, &IdealMqslsControl::poseInfoCallback, this)) {
        RCLCPP_ERROR(this->get_logger(), "Failed to subscribe to %s", pose_topic.c_str());
        return -1;
    }

    // publisher for EntityWrench msg
    std::string wrench_topic = "/world/" + world_name + "/virtual_wrench";
    _wrench_pub = _node.Advertise<gz::msgs::EntityWrench>(wrench_topic);

    RCLCPP_INFO(this->get_logger(), "Setup in %s", world_name.c_str());
    return 0;
}

void IdealMqslsControl::control_step(const uint64_t dt)
{
    // x500_1
    gz::msgs::EntityWrench wrench;
    wrench.mutable_entity()->set_name("x500_1");
    wrench.mutable_entity()->set_type(gz::msgs::Entity::MODEL);

    gz::math::Vector3d force{0, 0, 16.758};
    gz::math::Vector3d torque{0, 0, 0};

    gz::msgs::Set(wrench.mutable_wrench()->mutable_force(), force);
    gz::msgs::Set(wrench.mutable_wrench()->mutable_torque(), torque);

    // _wrench_pub.Publish(wrench);

    // // x500_2
    // wrench.mutable_entity()->set_name("x500_2");
    // _wrench_pub.Publish(wrench);

    // // x500_3
    // wrench.mutable_entity()->set_name("x500_3");
    // wrench.mutable_wrench()->mutable_force()->set_z(26.558);
    // _wrench_pub.Publish(wrench);

#if 1
    // simulate external disturbance
    wrench.mutable_entity()->set_name(_payload_name);
    wrench.mutable_entity()->set_type(gz::msgs::Entity::MODEL);
    if (_gz_sim_time > 5_s && _gz_sim_time < 5.5_s) {
        gz::msgs::Set(wrench.mutable_wrench()->mutable_force(), gz::math::Vector3d(0, 1, 0));
        gz::msgs::Set(wrench.mutable_wrench()->mutable_torque(), gz::math::Vector3d(0, 0, 0));
        RCLCPP_INFO(this->get_logger(), "External disturbance applied");
    }
    else {
        gz::msgs::Set(wrench.mutable_wrench()->mutable_force(), gz::math::Vector3d(0, 0, 0));
        gz::msgs::Set(wrench.mutable_wrench()->mutable_torque(), gz::math::Vector3d(0, 0, 0));
    }
    _wrench_pub.Publish(wrench);
#endif

    RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "Control force: %f", force.Z());
}

void IdealMqslsControl::clockCallback(const gz::msgs::Clock &clock)
{
    // get current time
    static uint64_t last_time = 0;
    _gz_sim_time = clock.sim().sec() * 1e6 + clock.sim().nsec() / 1e3;
    
    // dt [us] --> [1ms, 100ms]
    const uint64_t dt = std::clamp(_gz_sim_time - last_time, 1_ms, 100_ms);
    
    // control
    if (dt >= _control_interval) {
        // control logic
        control_step(dt);

        last_time = _gz_sim_time;
    }
}

void IdealMqslsControl::poseInfoCallback(const gz::msgs::Pose_V &pose)
{
#define get_ned_motion(ref_last_time, ref_postion, ref_velocity) \
{ \
    const uint64_t current_time = pose.header().stamp().sec() * 1e6 + pose.header().stamp().nsec() / 1e3; \
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
    ref_postion = position; \
    ref_velocity = velocity; \
}

    for (int i = 0; i < pose.pose_size(); i++) {
        if (pose.pose(i).name() == _payload_name) {
            
            static uint64_t last_time = 0;
            get_ned_motion(last_time, _load_position, _load_velocity);

            record(pose.header().stamp().sec() * 1e6 + pose.header().stamp().nsec() / 1e3);

            RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "%s Position: [%f, %f, %f], Velocity: [%f %f %f]", 
                            pose.pose(i).name().c_str(),
                            _load_position.x(), _load_position.y(), _load_position.z(),
                            _load_velocity.x(), _load_velocity.y(), _load_velocity.z());
        }
        else if (pose.pose(i).name().compare(0, _uav_name_prefix.size(), _uav_name_prefix) == 0) {
            int uav_id = std::stoi(pose.pose(i).name().substr(_uav_name_prefix.size())) - 1;
            if (uav_id < 0 || uav_id >= 3) {
                RCLCPP_ERROR(this->get_logger(), "Invalid UAV ID: %d", uav_id);
                return;
            }

            static uint64_t last_time[3] = {0};
            get_ned_motion(last_time[uav_id], _uav_position[uav_id], _uav_velocity[uav_id]);

            RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "%s Position: [%f, %f, %f], Velocity: [%f %f %f]", 
                            pose.pose(i).name().c_str(), 
                            _uav_position[uav_id].x(), _uav_position[uav_id].y(), _uav_position[uav_id].z(),
                            _uav_velocity[uav_id].x(), _uav_velocity[uav_id].y(), _uav_velocity[uav_id].z());
        }
    }

#undef get_ned_motion
}

} // namespace mqsls

int main(int argc, const char** argv) {
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<mqsls::IdealMqslsControl>(20_ms));

    rclcpp::shutdown();
    
    return 0;
}