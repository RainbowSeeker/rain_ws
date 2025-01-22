#include <ros/ros.h>
#include <Eigen/Eigen>
#include <geometry_msgs/PoseStamped.h>

#include <boost/asio.hpp>
#include <iostream>

struct MocapData
{
    int id; // 1, 2, 3 for uav; 4 for payload
    uint64_t timestamp;
    float position[3];  // Motion capture data
    float q[4];
};

using namespace boost::asio;

class Ros1MocapBridge
{
public:
    Ros1MocapBridge(const std::string& dst_ip, int port)
        : _remote_endpoint(ip::address::from_string(dst_ip), port),
          _socket(_io, ip::udp::endpoint(ip::udp::v4(), 0))
    {
        ROS_INFO("ROS1 Mocap Bridge starting");

        // Subscriber
        const std::string rigid_body_prefix = "rb_";
        for (int i = 0; i < 4; i++)
        {
            _odom_sub[i] = _nh.subscribe<geometry_msgs::PoseStamped>(
                "/vrpn_client_node/" + rigid_body_prefix + std::to_string(i + 1) + "/pose", 10, 
                [this, i](const geometry_msgs::PoseStamped::ConstPtr& msg) 
                {
                    pose_callback(msg, i + 1);
                });
        }
    }

    ~Ros1MocapBridge()
    {
    }

private:
    void rotateQuaternion(Eigen::Quaterniond &q_FRD_to_NED, const Eigen::Quaterniond q_FLU_to_ENU)
    {
        // FLU (ROS) to FRD (PX4) static rotation
        static const auto q_FLU_to_FRD = Eigen::Quaterniond(0, 1, 0, 0);

        /**
         * @brief Quaternion for rotation between ENU and NED frames
         *
         * NED to ENU: +PI/2 rotation about Z (Down) followed by a +PI rotation around X (old North/new East)
         * ENU to NED: +PI/2 rotation about Z (Up) followed by a +PI rotation about X (old East/new North)
         * This rotation is symmetric, so q_ENU_to_NED == q_NED_to_ENU.
         */
        static const auto q_ENU_to_NED = Eigen::Quaterniond(0, 0.70711, 0.70711, 0);

        // final rotation composition
        q_FRD_to_NED = q_ENU_to_NED * q_FLU_to_ENU * q_FLU_to_FRD.inverse();
    }

    void pose_callback(const geometry_msgs::PoseStamped::ConstPtr& msg, int id)
    {
        MocapData data;
        data.id = id;
        // data.timestamp = ros::Time::now().toNSec() / 1000;
        data.timestamp = msg->header.stamp.toNSec() / 1000;

        // ENU -> NED
        // TODO: Check the rotation
        data.position[0] = msg->pose.position.y;
        data.position[1] = msg->pose.position.x;
        data.position[2] = -msg->pose.position.z;
        Eigen::Quaterniond q_gr = Eigen::Quaterniond(
                                    msg->pose.orientation.w, 
                                    msg->pose.orientation.x,
                                    msg->pose.orientation.y,
                                    msg->pose.orientation.z);
        Eigen::Quaterniond q_nb;
        rotateQuaternion(q_nb, q_gr);

        data.q[0] = q_nb.w();
        data.q[1] = q_nb.x();
        data.q[2] = q_nb.y();
        data.q[3] = q_nb.z();

        _socket.send_to(buffer(&data, sizeof(data)), _remote_endpoint);

        static int count = 0;
        if (count ++ % 400 == 0)
            ROS_INFO("Send mocap data: id=%d, x=%f, y=%f, z=%f", data.id, data.position[0], data.position[1], data.position[2]);
    }

    // ros1
    ros::NodeHandle _nh;
    ros::Subscriber _odom_sub[4];

    // asio
    io_service _io;
    ip::udp::socket _socket;
    ip::udp::endpoint _remote_endpoint;
};

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: ros1_mocap_bridge <dst_ip> <port>" << std::endl;
        return 1;
    }

    ros::init(argc, argv, "ros1_mocap_bridge");

    Ros1MocapBridge bridge(argv[1], std::stoi(argv[2]));

    ros::spin();

    return 0;
}


/* TEST:
rostopic pub -r 100 /vrpn_client_node/rb_4/pose geometry_msgs/PoseStamped \
"header:
  seq: 0
  stamp:
    secs: 0
    nsecs: 0
  frame_id: 'world'
pose:
  position:
    x: 1.0
    y: 2.0
    z: 3.0
  orientation:
    x: 0.0
    y: 0.0
    z: 0.0
    w: 1.0"
*/

