#include <ros/ros.h>
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
    void pose_callback(const geometry_msgs::PoseStamped::ConstPtr& msg, int id)
    {
        MocapData data;
        data.id = id;
        // data.timestamp = ros::Time::now().toNSec() / 1000;
        data.timestamp = msg->header.stamp.toNSec() / 1000;
        data.position[0] = msg->pose.position.x;
        data.position[1] = msg->pose.position.y;
        data.position[2] = msg->pose.position.z;
        data.q[0] = msg->pose.orientation.w;
        data.q[1] = msg->pose.orientation.x;
        data.q[2] = msg->pose.orientation.y;
        data.q[3] = msg->pose.orientation.z;

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

