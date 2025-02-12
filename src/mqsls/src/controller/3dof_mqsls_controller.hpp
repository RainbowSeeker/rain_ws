#include <rclcpp/rclcpp.hpp>
#include <Eigen/Eigen>
#include <control_toolbox/pid.hpp>

#include <formation/data_recorder.hpp>
#include <formation/utils.hpp>
#include <mqsls/AlphaFilter.hpp>
#include <mqsls/px4_actuator.hpp>
#include <mqsls/trajectory_generator.hpp>
#include <mqsls/msg/follower_send.hpp>
#include <mqsls/msg/follower_recv.hpp>
#include <mqsls/srv/force_opt.hpp>
#include "model/controller.hpp"

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

    // ext state
    Eigen::Vector3d dL;
    double margin;

    Eigen::Vector3d load_position_sp;
    Eigen::Vector3d load_velocity_sp;

    // operator<<
    friend std::ostream &operator<<(std::ostream &os, const MqslsDataFrame &frame) {
        os << frame.timestamp << ' '
           << frame.load_position.x() << ' ' << frame.load_position.y() << ' ' << frame.load_position.z() << ' '
           << frame.load_velocity.x() << ' ' << frame.load_velocity.y() << ' ' << frame.load_velocity.z() << ' ';
        for (int i = 0; i < 3; i++) {
            os << frame.uav_position[i].x() << ' ' << frame.uav_position[i].y() << ' ' << frame.uav_position[i].z() << ' '
               << frame.uav_velocity[i].x() << ' ' << frame.uav_velocity[i].y() << ' ' << frame.uav_velocity[i].z() << ' ';
        }
        os << frame.dL.x() << ' ' << frame.dL.y() << ' ' << frame.dL.z() << ' ';
        os << frame.margin << ' ';
        os << frame.load_position_sp.x() << ' ' << frame.load_position_sp.y() << ' ' << frame.load_position_sp.z() << ' '
           << frame.load_velocity_sp.x() << ' ' << frame.load_velocity_sp.y() << ' ' << frame.load_velocity_sp.z() << ' ';
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
        ss  << "dLx dLy dLz" << ' '
            << "margin" << ' ';
        ss  << "x_sp y_sp z_sp" << ' '
            << "vx_sp vy_sp vz_sp" << ' ';

        return ss.str();
    }
};



} // namespace mqsls