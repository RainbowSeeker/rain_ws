// #include <rclcpp/rclcpp.hpp>
// #include <px4_msgs/msg/sensor_combined.hpp>

// /**
//  * @brief Sensor Combined uORB topic data callback
//  */
// class SensorCombinedListener : public rclcpp::Node
// {
// public:
// 	explicit SensorCombinedListener() : Node("sensor_combined_listener")
// 	{
// 		rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
// 		auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);
		
// 		subscription_ = this->create_subscription<px4_msgs::msg::SensorCombined>("/fmu/out/sensor_combined", qos,
// 		[this](const px4_msgs::msg::SensorCombined::UniquePtr msg) {
// 			std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
// 			std::cout << "RECEIVED SENSOR COMBINED DATA"   << std::endl;
// 			std::cout << "============================="   << std::endl;
// 			std::cout << "ts: "          << msg->timestamp    << std::endl;
// 			std::cout << "gyro_rad[0]: " << msg->gyro_rad[0]  << std::endl;
// 			std::cout << "gyro_rad[1]: " << msg->gyro_rad[1]  << std::endl;
// 			std::cout << "gyro_rad[2]: " << msg->gyro_rad[2]  << std::endl;
// 			std::cout << "gyro_integral_dt: " << msg->gyro_integral_dt << std::endl;
// 			std::cout << "accelerometer_timestamp_relative: " << msg->accelerometer_timestamp_relative << std::endl;
// 			std::cout << "accelerometer_m_s2[0]: " << msg->accelerometer_m_s2[0] << std::endl;
// 			std::cout << "accelerometer_m_s2[1]: " << msg->accelerometer_m_s2[1] << std::endl;
// 			std::cout << "accelerometer_m_s2[2]: " << msg->accelerometer_m_s2[2] << std::endl;
// 			std::cout << "accelerometer_integral_dt: " << msg->accelerometer_integral_dt << std::endl;
// 		});
// 	}

// private:
// 	rclcpp::Subscription<px4_msgs::msg::SensorCombined>::SharedPtr subscription_;
// };


// int main(int argc, char *argv[])
// {
// 	std::cout << "Starting sensor_combined listener node..." << std::endl;
// 	setvbuf(stdout, NULL, _IONBF, BUFSIZ);
// 	rclcpp::init(argc, argv);
// 	rclcpp::spin(std::make_shared<SensorCombinedListener>());

// 	rclcpp::shutdown();
// 	return 0;
// }

#include <chrono>
#include <rclcpp/rclcpp.hpp>
#include <px4_msgs/msg/debug_vect.hpp>

using namespace std::chrono_literals;

class DebugVectAdvertiser : public rclcpp::Node
{
public:
	DebugVectAdvertiser() : Node("debug_vect_advertiser") {
		publisher_ = this->create_publisher<px4_msgs::msg::DebugVect>("fmu/debug_vect/in", 10);
		auto timer_callback =
		[this]()->void {
			auto debug_vect = px4_msgs::msg::DebugVect();
			debug_vect.timestamp = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::steady_clock::now()).time_since_epoch().count();
			std::string name = "test";
			std::copy(name.begin(), name.end(), debug_vect.name.begin());
			debug_vect.x = 1.0;
			debug_vect.y = 2.0;
			debug_vect.z = 3.0;
			RCLCPP_INFO(this->get_logger(), "\033[97m Publishing debug_vect: time: %llu x: %f y: %f z: %f \033[0m",
                                debug_vect.timestamp, debug_vect.x, debug_vect.y, debug_vect.z);
			this->publisher_->publish(debug_vect);
		};
		timer_ = this->create_wall_timer(500ms, timer_callback);

		std::cout << "clk " << now() << std::endl;
	}

private:
	rclcpp::TimerBase::SharedPtr timer_;
	rclcpp::Publisher<px4_msgs::msg::DebugVect>::SharedPtr publisher_;
};

int main(int argc, char *argv[])
{
	std::cout << "Starting debug_vect advertiser node..." << std::endl;
	setvbuf(stdout, NULL, _IONBF, BUFSIZ);
	rclcpp::init(argc, argv);
	rclcpp::spin(std::make_shared<DebugVectAdvertiser>());

	rclcpp::shutdown();
	return 0;
}