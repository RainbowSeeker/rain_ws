#include <chrono>
#include <rclcpp/rclcpp.hpp>
#include <px4_msgs/msg/vehicle_local_position.hpp>
#include <px4_msgs/msg/formation_cross.hpp>

using namespace std::chrono_literals;

class cross_node : public rclcpp::Node
{
public:
    cross_node() : Node("cross_node")
    {
        rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
		auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);
		
        for (int i = 0; i < 3; i++)
        {
            // Subscribe to the vehicle_local_position topic
            _local_position_sub[i] = this->create_subscription<px4_msgs::msg::VehicleLocalPosition>(
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

            _formation_cross_pub[i] = this->create_publisher<px4_msgs::msg::FormationCross>(
            std::string("px4_") + (char)('1' + i) + "/fmu/in/formation_cross", 10);
        }
        
        auto timer_callback = [this]() -> void {
            for (int i = 0; i < 3; i++)
                _formation_cross_pub[i]->publish(_formation_cross);
        };

        _timer = this->create_wall_timer(40ms, timer_callback);
    }
    ~cross_node() = default;

private:
    rclcpp::TimerBase::SharedPtr _timer;
    rclcpp::Subscription<px4_msgs::msg::VehicleLocalPosition>::SharedPtr _local_position_sub[3];
    rclcpp::Publisher<px4_msgs::msg::FormationCross>::SharedPtr _formation_cross_pub[3];

    px4_msgs::msg::FormationCross _formation_cross{};
};

int main(int argc, char *argv[])
{
	rclcpp::init(argc, argv);
	rclcpp::spin(std::make_shared<cross_node>());

	rclcpp::shutdown();
	return 0;
}