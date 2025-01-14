#ifndef GZ_COMPONENT_HPP
#define GZ_COMPONENT_HPP

#include "component.hpp"
#include "gz/msgs.hh"
#include "gz/transport.hh"
#include "gz/math.hh"
#include <iostream>
#include <string>

namespace mqsls {

class GzInputSource : public InputSource
{
public:
    const uint32_t rope_id_map[3] = {49, 172, 221};
    const uint32_t rope_id_offset = 7;

    GzInputSource(gz::transport::Node *node, int node_index, const std::string &world_name)
     : _node(node), _uav_name("x500_" + std::to_string(node_index)), _rope_id(rope_id_map[node_index - 1])
    {
        // pose: /world/$WORLD/pose/info
        std::string pose_topic = "/world/" + world_name + "/pose/info";
        if (!_node->Subscribe(pose_topic, &GzInputSource::poseInfoCallback, this)) {
            throw std::runtime_error("Failed to subscribe to " + pose_topic);
        }

        // IMU: /world/$WORLD/model/rope_xxx/link/link_5/sensor/imu_sensor/imu
        std::string imu_topic = "/world/" + world_name + "/model/rope_" + std::to_string(node_index) + "/link/link_5/sensor/imu_sensor/imu";
        if (!_node->Subscribe(imu_topic, &GzInputSource::imuCallback, this)) {
            throw std::runtime_error("Failed to subscribe to " + imu_topic);
        }
    }

    void register_update_callback(std::function<void(const InputData &)> callback) override
    {
        _callbacks.push_back(callback);
    }
    
private:
    const InputData &convert(uint64_t timestamp)
    {
        _data.msg.timestamp = timestamp;
        _data.msg.position_uav[0] = _uav_position.x();
        _data.msg.position_uav[1] = _uav_position.y();
        _data.msg.position_uav[2] = _uav_position.z();
        _data.msg.velocity_uav[0] = _uav_velocity.x();
        _data.msg.velocity_uav[1] = _uav_velocity.y();
        _data.msg.velocity_uav[2] = _uav_velocity.z();

        // _data.msg.delta_position[0] = _cable_direction[0] * 2;
        // _data.msg.delta_position[1] = _cable_direction[1] * 2;
        // _data.msg.delta_position[2] = _cable_direction[2] * 2;
        _data.msg.position_load[0] = _load_position.x();
        _data.msg.position_load[1] = _load_position.y();
        _data.msg.position_load[2] = _load_position.z();
        _data.msg.velocity_load[0] = _load_velocity.x();
        _data.msg.velocity_load[1] = _load_velocity.y();
        _data.msg.velocity_load[2] = _load_velocity.z();

        return _data;
    }

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
                
                static uint64_t last_time = 0;
                get_ned_motion(last_time, _load_position, _load_velocity);

                // saturation
                _load_velocity = _load_velocity.cwiseMax(min_velocity).cwiseMin(max_velocity);

            }
            else if (pose.pose(i).name() == _uav_name) {

                static uint64_t last_time = {0};
                get_ned_motion(last_time, _uav_position, _uav_velocity);

                // saturation
                _uav_velocity= _uav_velocity.cwiseMax(min_velocity).cwiseMin(max_velocity);

                for (const auto &cb : _callbacks) {
                    cb(convert(current_time));
                }
            }
            else if (pose.pose(i).id() == _rope_id) {
                _cable_att_enu = {
                    pose.pose(i).orientation().w(),
                    pose.pose(i).orientation().x(),
                    pose.pose(i).orientation().y(),
                    pose.pose(i).orientation().z()
                };
            }
            else if (pose.pose(i).id() == _rope_id + rope_id_offset) {
                static uint64_t last_time = 0;
                const double dt = std::clamp((current_time - last_time) * 1e-6, 0.001, 0.1);
                last_time = current_time;

                Eigen::Quaterniond relative_att = {
                    pose.pose(i).orientation().w(),
                    pose.pose(i).orientation().x(),
                    pose.pose(i).orientation().y(),
                    pose.pose(i).orientation().z()
                };

                const Eigen::Vector3d flu_direction = _cable_att_enu * relative_att * Eigen::Vector3d(0, 0, -1);
                const Eigen::Vector3d frd_direction = {flu_direction.y(), flu_direction.x(), -flu_direction.z()};
                const Eigen::Vector3d dot_direction = (frd_direction - _cable_direction) / dt;
                
                _cable_direction = frd_direction;

                static const Eigen::Vector3d max_dot_direction{0.5, 0.5, 0.5}; 
                static const Eigen::Vector3d min_dot_direction{-0.5, -0.5, -0.5};
                _cable_dot_direction = dot_direction.cwiseMax(min_dot_direction).cwiseMin(max_dot_direction);
            }
        }

    #undef get_ned_motion
    }

    void imuCallback(const gz::msgs::IMU &imu)
    {
        // FLU -> FRD
        static const auto q_FLU_to_FRD = gz::math::Quaterniond(0, 1, 0, 0);

        gz::math::Vector3d accel_b = q_FLU_to_FRD.RotateVector(gz::math::Vector3d(
					     imu.linear_acceleration().x(),
					     imu.linear_acceleration().y(),
					     imu.linear_acceleration().z()));

        gz::math::Vector3d gyro_b = q_FLU_to_FRD.RotateVector(gz::math::Vector3d(
					    imu.angular_velocity().x(),
					    imu.angular_velocity().y(),
					    imu.angular_velocity().z()));
    }

    InputData _data;
    gz::transport::Node *_node;
    const std::string _payload_name = "payload";
    const std::string _uav_name;
    const uint32_t _rope_id;
    std::vector<std::function<void(const InputData &)>> _callbacks;
    Eigen::Vector3d _load_position;
    Eigen::Vector3d _load_velocity;
    Eigen::Vector3d _uav_position;
    Eigen::Vector3d _uav_velocity;

    // addon
    Eigen::Quaterniond _cable_att_enu;
    Eigen::Vector3d _cable_direction {};
    Eigen::Vector3d _cable_dot_direction {};
};

class GzOutputActuator : public OutputActuator
{
public:
    GzOutputActuator(gz::transport::Node *node, int node_index, const std::string &world_name) : _node(node), _uav_name("x500_" + std::to_string(node_index))
    {
        // publisher for EntityWrench msg
        std::string wrench_topic = "/world/" + world_name + "/virtual_wrench";
        _wrench_pub = _node->Advertise<gz::msgs::EntityWrench>(wrench_topic);
    }

    void apply(const OutputData &data) override
    {
        // x500_i
        gz::msgs::EntityWrench wrench;
        wrench.mutable_entity()->set_name(_uav_name);
        wrench.mutable_entity()->set_type(gz::msgs::Entity::MODEL);
        wrench.mutable_wrench()->mutable_force()->set_x(data.msg.force[1]);
        wrench.mutable_wrench()->mutable_force()->set_y(data.msg.force[0]);
        wrench.mutable_wrench()->mutable_force()->set_z(-data.msg.force[2]);

        _wrench_pub.Publish(wrench);
    }

private:
    gz::transport::Node *_node;
    const std::string _uav_name;
    gz::transport::Node::Publisher _wrench_pub;
};

class GzEventHandler : public EventHandler
{
public:
    GzEventHandler(gz::transport::Node *node, const std::string &world_name) : _node(node)
    {    
        // clock
        std::string clock_topic = "/world/" + world_name + "/clock";
        if (!_node->Subscribe(clock_topic, &GzEventHandler::clockCallback, this)) {
            std::cerr << "Failed to subscribe to " << clock_topic << std::endl;
            return;
        }

        // stats
        std::string stats_topic = "/world/" + world_name + "/stats";
        if (!_node->Subscribe(stats_topic, &GzEventHandler::statsCallback, this)) {
            std::cerr << "Failed to subscribe to " << stats_topic << std::endl;
            return;
        }
        
    }

    void register_periodic_callback(std::function<void(uint64_t)> callback, uint64_t period_us) override
    {
        _periodic_callbacks.push_back({callback, period_us, 0});
    }

private:
    // callbacks
    void statsCallback(const gz::msgs::WorldStatistics &stats)
    {
        static bool first = true;
        bool paused = stats.paused();
        if (first || paused != _is_paused) {
            std::cout << "Simulation is " << (paused ? "paused" : "running") << std::endl;
            first = false;
        }
        _is_paused = paused;
    }

    void clockCallback(const gz::msgs::Clock &clock)
    {
        if (_is_paused) {
            return;
        }
        // get current time
        _gz_sim_time = clock.sim().sec() * 1e6 + clock.sim().nsec() / 1e3;

        // loop
        for (auto &cb : _periodic_callbacks) {
            if (_gz_sim_time - cb.last_time >= cb.period_us) {
                cb.callback(_gz_sim_time);
                cb.last_time = _gz_sim_time;
            }
        }
    }

    gz::transport::Node *_node;
    struct cb_body
    {
        std::function<void(uint64_t)> callback;
        uint64_t period_us;
        uint64_t last_time;
    };
    std::vector<cb_body> _periodic_callbacks; 
    uint64_t _gz_sim_time {0}; // [us]
    bool    _is_paused {true};
};

} // namespace mqsls

#endif // !GZ_COMPONENT_HPP