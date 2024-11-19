#include <rclcpp/rclcpp.hpp>
#include <px4_msgs/msg/vehicle_command.hpp>

int main(int argc, const char** argv) {

    if (argc < 3) {
        std::cerr << "Usage: commander <obj_id> <command> [<param1> <param2> ...]" << std::endl;
        return 1;
    }

    int obj_id = std::stoi(argv[1]);
    uint32_t command = static_cast<uint32_t>(std::stoul(argv[2]));
    float param1 = argc > 3 ? std::stof(argv[3]) : NAN;
    float param2 = argc > 4 ? std::stof(argv[4]) : NAN;
    float param3 = argc > 5 ? std::stof(argv[5]) : NAN;
    float param4 = argc > 6 ? std::stof(argv[6]) : NAN;
    double param5 = argc > 7 ? std::stod(argv[7]) : NAN;
    double param6 = argc > 8 ? std::stod(argv[8]) : NAN;
    float param7 = argc > 9 ? std::stof(argv[9]) : NAN;

    rclcpp::init(argc, argv);

    std::string topic_ns = obj_id > 0 ? "/px4_" + std::to_string(obj_id) : "";
    auto node = rclcpp::Node::make_shared("commander");
    auto publisher = node->create_publisher<px4_msgs::msg::VehicleCommand>(topic_ns + "/fmu/in/vehicle_command", 10);
    auto cmd = px4_msgs::msg::VehicleCommand();
    cmd.command = command;
    cmd.param1 = param1;
    cmd.param2 = param2;
    cmd.param3 = param3;
    cmd.param4 = param4;
    cmd.param5 = param5;
    cmd.param6 = param6;
    cmd.param7 = param7;
    cmd.target_system = obj_id;
    cmd.source_system = obj_id;
    cmd.target_component = 0;
    cmd.source_component = 0;
    cmd.from_external = true;
    cmd.timestamp = node->get_clock()->now().nanoseconds() / 1e3;
    publisher->publish(cmd);

    rclcpp::shutdown();

    return 0;
}