#pragma once
#include <rclcpp/rclcpp.hpp>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include "component.hpp"

namespace mqsls {

struct MqslsDataFrame {
    // timestamp
    uint64_t timestamp;
    // payload state
    Eigen::Vector3d load_position;
    Eigen::Vector3d load_velocity;
    Eigen::Vector3d load_euler;
    Eigen::Vector3d load_omega;

    // cable state
    Eigen::Vector3d q[3];
    Eigen::Vector3d w[3];

    // operator<<
    friend std::ostream &operator<<(std::ostream &os, const MqslsDataFrame &frame) {
        os << frame.timestamp << ' '
           << frame.load_position.x() << ' ' << frame.load_position.y() << ' ' << frame.load_position.z() << ' '
           << frame.load_velocity.x() << ' ' << frame.load_velocity.y() << ' ' << frame.load_velocity.z() << ' '
           << frame.load_euler.x() << ' ' << frame.load_euler.y() << ' ' << frame.load_euler.z() << ' '
           << frame.load_omega.x() << ' ' << frame.load_omega.y() << ' ' << frame.load_omega.z() << ' ';
        for (int i = 0; i < 3; i++) {
            os << frame.q[i].x() << ' ' << frame.q[i].y() << ' ' << frame.q[i].z() << ' '
               << frame.w[i].x() << ' ' << frame.w[i].y() << ' ' << frame.w[i].z() << ' ';
        }
        os << '\n';
        return os;
    }

    static std::string header() {
        std::stringstream ss;
        ss  << "timestamp" << ' '
            << "x y z" << ' '
            << "vx vy vz" << ' '
            << "phi theta psi" << ' '
            << "p q r" << ' ';
        for (int i = 1; i <= 3; i++) {
            ss  << "qx" << i << " qy" << i << " qz" << i << ' '
                << "wx" << i << " wy" << i << " wz" << i << ' ';
        }

        return ss.str();
    }
};


} // namespace mqsls