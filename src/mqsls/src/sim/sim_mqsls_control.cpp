#include <control_toolbox/pid.hpp>

#include "formation/data_recorder.hpp"
#include "formation/utils.hpp"
#include "../model/interface.hpp"
#include "mqsls/mqsls.hpp"

using namespace std::chrono_literals;

namespace mqsls {

static gz::transport::Node &get_default_gz_node()
{
    static gz::transport::Node node;
    return node;
}

static std::shared_ptr<InputSource> make_input_source(rclcpp::Node *node, int node_index, const std::string &input_source, const std::string &world_name)
{
    if (input_source == "gz") {
        return std::make_shared<GzInputSource>(&get_default_gz_node(), node_index, world_name);
    } else if (input_source == "px4") {
        return std::make_shared<PX4InputSource>();
    } else {
        throw std::runtime_error("Invalid input source");
    }
}

static std::shared_ptr<OutputActuator> make_output_actuator(rclcpp::Node *node, int node_index, const std::string &output_actuator, const std::string &world_name)
{
    if (output_actuator == "gz") {
        return std::make_shared<GzOutputActuator>(&get_default_gz_node(), node_index, world_name);
    } else if (output_actuator == "px4") {
        return std::make_shared<PX4OutputActuator>(node, node_index);
    } else {
        throw std::runtime_error("Invalid output actuator");
    }
}

static std::shared_ptr<EventHandler> make_event_handler(rclcpp::Node *node, const std::string &event_handler, const std::string &world_name)
{
    if (event_handler == "gz") {
        return std::make_shared<GzEventHandler>(&get_default_gz_node(), world_name);
    } else if (event_handler == "px4") {
        return std::make_shared<PX4EventHandler>(node);
    } else {
        throw std::runtime_error("Invalid event handler");
    }
}

static void delta_pv_to_qw(const Eigen::Vector3d &delta_position, const Eigen::Vector3d &delta_velocity, Eigen::Vector3d &q, Eigen::Vector3d &w)
{
    q = delta_position.normalized();
    w = q.cross(delta_velocity / delta_position.norm() - delta_position * delta_position.dot(delta_velocity) / std::pow(delta_position.norm(), 3));
}

static void cable_direction_planner(Eigen::Vector3d q[3])
{
    #define deg2rad(x) ((x) * M_PI / 180)
    const Eigen::Vector3d theta_sp = {deg2rad(60), deg2rad(60), deg2rad(60)};
    const Eigen::Vector3d psi_sp = {deg2rad(-120), deg2rad(0), deg2rad(120)};

    for (int i = 0; i < 3; i++) {
        q[i] = Eigen::Vector3d(-cos(psi_sp[i]) * cos(theta_sp[i]), -sin(psi_sp[i]) * cos(theta_sp[i]), sin(theta_sp[i]));
    }
}

class MqslsLeader : public rclcpp::Node
{
public:
    MqslsLeader(int node_index) : rclcpp::Node("amc_" + std::to_string(node_index))
    {
        parameter_declare();

        const std::string world_name = this->get_parameter("world_name").as_string();
        const std::string input_source = this->get_parameter("input_source").as_string();
        const std::string output_actuator = this->get_parameter("output_actuator").as_string();
        const std::string event_handler = this->get_parameter("event_handler").as_string();
        const double boot_wait = this->get_parameter("boot_wait").as_double();

        RCLCPP_INFO(this->get_logger(), "Node %d boot wait: %f", node_index, boot_wait);
        std::this_thread::sleep_for(std::chrono::microseconds(static_cast<int>(boot_wait * 1e6)));
        
        _input_source = make_input_source(this, node_index, input_source, world_name);
        _output_actuator = make_output_actuator(this, node_index, output_actuator, world_name); // zero-hold actuator
        _event_handler = make_event_handler(this, event_handler, world_name);

        rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
        auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);

        for (int i = 0; i < 2; i++) {
            _leader_recv_sub[i] = this->create_subscription<mqsls::msg::FollowerSend>(
                "follower_send" + std::to_string(i + 2), qos,
                [this, i](const mqsls::msg::FollowerSend::SharedPtr msg) {
                    _follower_msg[i + 1] = *msg;
                });
        }

        for (int i = 0; i < 2; i++) {
            _leader_send_pub[i] = this->create_publisher<mqsls::msg::FollowerRecv>("follower_recv" + std::to_string(i + 2), 10);
        }

        _input_source->register_update_callback([this](const InputSource::InputData &data) {
            _follower_msg[0] = data.msg;
        });

        _event_handler->register_periodic_callback(std::bind(&MqslsLeader::run, this, std::placeholders::_1), 20_ms);

        cable_direction_planner(_cable_dir_sp);
    }
    ~MqslsLeader() = default;
private:
    enum MQSLS_RUNNING_STATE
    {
        STATE_PREPARE,
        STATE_RUNNING,
        STATE_STOPPING,
        STATE_OVER
    } _state = STATE_PREPARE, _last_state = _state;

    void state_machine_flow(CodeGenController::OutputBus &output)
    {
        switch (_state)
        {
        case STATE_PREPARE:
            if (state_prepare(output)) {
                _state = STATE_RUNNING;
            }
            break;
        case STATE_RUNNING:
            state_running(output);
            break;
        case STATE_STOPPING:
            state_stopping(output);
            _state = STATE_OVER;
            break;
        case STATE_OVER:
            state_over(output);
            break;
        default:
            RCLCPP_ERROR(this->get_logger(), "Invalid state");
            break;
        }
        if (_state != _last_state) {
            RCLCPP_INFO(this->get_logger(), "State transition: %d -> %d", _last_state, _state);
            _last_state = _state;
        }
    }

    bool state_prepare(CodeGenController::OutputBus &output)
    {
        // prepare phase: set uav pose to expected cope end
        const double uav_mass = this->get_parameter("uav_mass").as_double();
        const double cable_len = this->get_parameter("cable_len").as_double();
        const double KP = 0.5;
        const double KV = 1;
        
        auto cable_end_controller = [&, this](int index, Eigen::Vector3d &position_err, Eigen::Vector3d &acceleration_sp) -> void
        {
            // NED frame, payload position is the origin
            const Eigen::Vector3d q_sp = _cable_dir_sp[index];
            const double safety_len = cable_len * 0.8; // 80% of cable length
            const Eigen::Vector3d position_sp = -q_sp * safety_len;
            const Eigen::Vector3d position_now = {-_follower_msg[index].delta_position[0], -_follower_msg[index].delta_position[1], -_follower_msg[index].delta_position[2]};
            const Eigen::Vector3d velocity_now = {_follower_msg[index].velocity_uav[0], _follower_msg[index].velocity_uav[1], _follower_msg[index].velocity_uav[2]};

            position_err = position_sp - position_now;
            acceleration_sp = KV * (KP * position_err - velocity_now);

            RCLCPP_DEBUG(this->get_logger(), "UAV%d perr: %f %f %f", index, position_err[0], position_err[1], position_err[2]);
        };

        Eigen::Vector3d acceleration_sp, position_err[3];
        #define FILL_OUTPUT_FORCE(i) \
        { \
            cable_end_controller(i - 1, position_err[i - 1], acceleration_sp); \
            output.force_sp##i[0] = uav_mass * acceleration_sp[0]; \
            output.force_sp##i[1] = uav_mass * acceleration_sp[1]; \
            output.force_sp##i[2] = uav_mass * acceleration_sp[2] - uav_mass * 9.81; \
        }
        
        FILL_OUTPUT_FORCE(1); // leader
        FILL_OUTPUT_FORCE(2); // follower 1
        FILL_OUTPUT_FORCE(3); // follower 2

        #undef FILL_OUTPUT_FORCE
    
        // check if all uav are in position
        for (int i = 0; i < 3; i++) {
            if (position_err[i].norm() > 0.3 || std::abs(position_err[i][2]) > 0.2) {
                return false;
            }
        }
        return true;
    }

    void state_running(CodeGenController::OutputBus &output)
    {
        // build input
        Eigen::Vector3d q[3], w[3];
        for (int i = 0; i < 3; i++) {
            delta_pv_to_qw(Eigen::Vector3d(_follower_msg[i].delta_position[0], _follower_msg[i].delta_position[1], _follower_msg[i].delta_position[2]),
                           Eigen::Vector3d(_follower_msg[i].delta_velocity[0], _follower_msg[i].delta_velocity[1], _follower_msg[i].delta_velocity[2]),
                           q[i], w[i]);
        }
        for (int i = 0; i < 3; i++) {
            _control_input.Payload_Out.pL[i] = _follower_msg[0].position_uav[i] + _follower_msg[0].delta_position[i];
            _control_input.Payload_Out.vL[i] = _follower_msg[0].velocity_uav[i] + _follower_msg[0].delta_velocity[i];
            _control_input.Payload_Out.q_1[i] = q[0][i];
            _control_input.Payload_Out.q_2[i] = q[1][i];
            _control_input.Payload_Out.q_3[i] = q[2][i];
            _control_input.Payload_Out.w_1[i] = w[0][i];
            _control_input.Payload_Out.w_2[i] = w[1][i];
            _control_input.Payload_Out.w_3[i] = w[2][i];
        }
        _control_input.Payload_Out.timestamp = _follower_msg[0].timestamp;
        
        // running
        _controller.step(_control_input);

        output = _controller.getOutput();
    }

    void state_stopping(CodeGenController::OutputBus &output)
    {
        // stopping
    }

    void state_over(CodeGenController::OutputBus &output)
    {
        // over
    }

    void run(uint64_t timestamp)
    {
        static uint64_t boot_time = timestamp;
        _running_time = timestamp - boot_time;

        parameter_update();

        CodeGenController::OutputBus output;

        state_machine_flow(output);
        
        // publish leader id == 1
        auto data = OutputActuator::OutputData();
        data.msg.timestamp = _running_time;
        data.msg.force[0] = output.force_sp1[0];
        data.msg.force[1] = output.force_sp1[1];
        data.msg.force[2] = output.force_sp1[2];
        _output_actuator->apply(data); // zero-hold actuator

        // publish follower id == 2
        mqsls::msg::FollowerRecv msg;
        msg.timestamp = _running_time;
        msg.force[0] = output.force_sp2[0];
        msg.force[1] = output.force_sp2[1];
        msg.force[2] = output.force_sp2[2];
        _leader_send_pub[0]->publish(msg);

        // publish follower id == 3
        msg.force[0] = output.force_sp3[0];
        msg.force[1] = output.force_sp3[1];
        msg.force[2] = output.force_sp3[2];
        _leader_send_pub[1]->publish(msg);

        record();
        
        RCLCPP_INFO(this->get_logger(), "\n UAV1: %f %f %f\n UAV2: %f %f %f\n UAV3: %f %f %f", 
                    output.force_sp1[0], output.force_sp1[1], output.force_sp1[2],
                    output.force_sp2[0], output.force_sp2[1], output.force_sp2[2],
                    output.force_sp3[0], output.force_sp3[1], output.force_sp3[2]);
    }

    // parameters
    void parameter_declare()
    {
        this->declare_parameter("world_name", "mqsls");
        this->declare_parameter("input_source", "gz"); // gz or px4
        this->declare_parameter("output_actuator", "gz"); // gz or px4
        this->declare_parameter("event_handler", "gz"); // gz or px4
        this->declare_parameter("boot_wait", 0.0); // boot wait time

        // controller parameters
        this->declare_parameter("cable_len", 1.0);
        this->declare_parameter("load_mass", 1.0);
        this->declare_parameter("uav_mass", 1.5);
    }

    void parameter_update()
    {
        CONTROL_PARAM.CABLE_LEN = this->get_parameter("cable_len").as_double();
        CONTROL_PARAM.MASS_LOAD = this->get_parameter("load_mass").as_double();
        CONTROL_PARAM.MASS_UAV = this->get_parameter("uav_mass").as_double();

        // CONTROL_PARAM.KP = 0.05;
        // CONTROL_PARAM.KV = 0.01;
        // CONTROL_PARAM.KQ = 1;
        // CONTROL_PARAM.KW = 2;
    }
    
    // Subscriber
    rclcpp::Subscription<mqsls::msg::FollowerSend>::SharedPtr _leader_recv_sub[2];


    // Publisher
    rclcpp::Publisher<mqsls::msg::FollowerRecv>::SharedPtr _leader_send_pub[2];

    // detail
    uint64_t _running_time = 0;
    Eigen::Vector3d _cable_dir_sp[3];
    mqsls::msg::FollowerSend _follower_msg[3]; // 0: leader, 1: follower 2, 2: follower 3
    std::shared_ptr<InputSource> _input_source;
    std::shared_ptr<OutputActuator> _output_actuator;
    std::shared_ptr<EventHandler> _event_handler;

    // controller
    CodeGenController::InputBus _control_input;
    CodeGenController _controller;

    // recorder
    void record()
    {
        MqslsDataFrame frame;
        frame.timestamp = _control_input.Payload_Out.timestamp;
        frame.load_position = Eigen::Vector3d(_control_input.Payload_Out.pL[0], _control_input.Payload_Out.pL[1], _control_input.Payload_Out.pL[2]);
        frame.load_velocity = Eigen::Vector3d(_control_input.Payload_Out.vL[0], _control_input.Payload_Out.vL[1], _control_input.Payload_Out.vL[2]);
        frame.q[0] = Eigen::Vector3d(_control_input.Payload_Out.q_1[0], _control_input.Payload_Out.q_1[1], _control_input.Payload_Out.q_1[2]);
        frame.q[1] = Eigen::Vector3d(_control_input.Payload_Out.q_2[0], _control_input.Payload_Out.q_2[1], _control_input.Payload_Out.q_2[2]);
        frame.q[2] = Eigen::Vector3d(_control_input.Payload_Out.q_3[0], _control_input.Payload_Out.q_3[1], _control_input.Payload_Out.q_3[2]);
        frame.w[0] = Eigen::Vector3d(_control_input.Payload_Out.w_1[0], _control_input.Payload_Out.w_1[1], _control_input.Payload_Out.w_1[2]);
        frame.w[1] = Eigen::Vector3d(_control_input.Payload_Out.w_2[0], _control_input.Payload_Out.w_2[1], _control_input.Payload_Out.w_2[2]);
        frame.w[2] = Eigen::Vector3d(_control_input.Payload_Out.w_3[0], _control_input.Payload_Out.w_3[1], _control_input.Payload_Out.w_3[2]);
        _recorder.push(frame);
    }
    DataRecorder<MqslsDataFrame> _recorder{"install/" + utils::nowstr() + ".csv", 10};
};

class MqslsFollower : public rclcpp::Node
{
public:
    MqslsFollower(int node_index) : rclcpp::Node("amc_" + std::to_string(node_index))
    {
        parameter_declare();

        const std::string world_name = this->get_parameter("world_name").as_string();
        const std::string input_source = this->get_parameter("input_source").as_string();
        const std::string output_actuator = this->get_parameter("output_actuator").as_string();
        _input_source = make_input_source(this, node_index, input_source, world_name);
        _output_actuator = make_output_actuator(this, node_index, output_actuator, world_name); // zero-hold actuator

        rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
        auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);

        _follower_recv_sub = this->create_subscription<mqsls::msg::FollowerRecv>(
            "follower_recv" + std::to_string(node_index), qos,
            [this](const mqsls::msg::FollowerRecv::SharedPtr msg) {
                auto data = OutputActuator::OutputData();
                data.msg = *msg;
                _output_actuator->apply(data); // zero-hold actuator
                RCLCPP_DEBUG(this->get_logger(), "Received force: %f %f %f", msg->force[0], msg->force[1], msg->force[2]);
            });
        
        _follower_send_pub = this->create_publisher<mqsls::msg::FollowerSend>("follower_send" + std::to_string(node_index), 10);

        _input_source->register_update_callback([this](const InputSource::InputData &data) {            
            _follower_send_pub->publish(data.msg);
        });
    }
    ~MqslsFollower() = default;
private:
    // parameters
    void parameter_declare()
    {
        this->declare_parameter("world_name", "mqsls");
        this->declare_parameter("input_source", "gz"); // gz or px4
        this->declare_parameter("output_actuator", "gz"); // gz or px4
    }
    // Subscriber
    rclcpp::Subscription<mqsls::msg::FollowerRecv>::SharedPtr _follower_recv_sub;

    // Publisher
    rclcpp::Publisher<mqsls::msg::FollowerSend>::SharedPtr _follower_send_pub;

    // detail
    std::shared_ptr<InputSource> _input_source;
    std::shared_ptr<OutputActuator> _output_actuator;
};
} // namespace mqsls


int main(int argc, const char** argv) 
{
    rclcpp::init(argc, argv);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <node_index>" << std::endl;
        return -1;
    }
    
    int index = std::stoi(argv[1]);
    // choose leader or follower instance to run
    switch (index)
    {
    case 1:
        rclcpp::spin(std::make_shared<mqsls::MqslsLeader>(index));
        break;
    case 2:
    case 3:
        rclcpp::spin(std::make_shared<mqsls::MqslsFollower>(index));
        break;
    default:
        throw std::runtime_error("Invalid node index");
    }

    rclcpp::shutdown();
    
    return 0;
}