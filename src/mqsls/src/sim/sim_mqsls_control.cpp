#include <control_toolbox/pid.hpp>

#include "formation/data_recorder.hpp"
#include "formation/utils.hpp"
#include "../model/interface.hpp"
#include "mqsls/mqsls.hpp"
#include "mqsls/AlphaFilter.hpp"
#include "mqsls/trajectory_generator.hpp"

#include <mqsls/srv/force_opt.hpp>

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

#define deg2rad(x) ((x) * M_PI / 180)

static std::shared_ptr<TrajectoryGenerator> make_trajectory_generator(const std::string &traj_type)
{
    if (traj_type == "line") {
        return std::make_shared<LineTrajectoryGenerator>(Eigen::Vector3d(0, 0, -10), Eigen::Vector3d(100, 0, -10), 2.0);
    } else if (traj_type == "circle") {
        return std::make_shared<CircleTrajectoryGenerator>(Eigen::Vector3d(0, 0, -10), 10, deg2rad(6));
    } else if (traj_type == "rectangle") {
        return std::make_shared<RectangleTrajectoryGenerator>(Eigen::Vector3d(0, 0, -10), Eigen::Vector3d(20, 20, -10), 2);
    } else {
        throw std::runtime_error("Invalid trajectory type");
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
        const std::string traj_type = this->get_parameter("traj_type").as_string();
        const double boot_wait = this->get_parameter("boot_wait").as_double();

        RCLCPP_INFO(this->get_logger(), "Node %d boot wait: %f", node_index, boot_wait);
        std::this_thread::sleep_for(std::chrono::microseconds(static_cast<int>(boot_wait * 1e6)));
        
        _input_source = make_input_source(this, node_index, input_source, world_name);
        _output_actuator = make_output_actuator(this, node_index, output_actuator, world_name); // zero-hold actuator
        _event_handler = make_event_handler(this, event_handler, world_name);
        _traj_gen = make_trajectory_generator(traj_type);

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

        _force_opt_client = this->create_client<mqsls::srv::ForceOpt>("force_opt");
        while (!_force_opt_client->wait_for_service(1s)) {
            if (!rclcpp::ok()) {
                RCLCPP_ERROR(this->get_logger(), "Interrupted while waiting for the service. Exiting.");
                return;
            }
            RCLCPP_INFO(this->get_logger(), "service not available, waiting again...");
        }
        
        RCLCPP_INFO(this->get_logger(), "Requesting initial cable direction");

        async_request_force_opt({0, 0, -_load_mass * 9.8});

        // start
        _event_handler->register_periodic_callback(std::bind(&MqslsLeader::run, this, std::placeholders::_1), 20_ms);
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
        const double KP = 1;
        const double KD = 1.2;
        
        auto cable_end_controller = [&, this](int index, Eigen::Vector3d &position_err, Eigen::Vector3d &acceleration_sp) -> void
        {
            // NED frame, payload position is the origin
            const Eigen::Vector3d q_sp = _cable_dir_sp[index];
            const double safety_len = _cable_len * 0.9; // 90% of cable length
            const Eigen::Vector3d position_sp = -q_sp * safety_len;
            const Eigen::Vector3d position_now = {-_follower_msg[index].delta_position[0], -_follower_msg[index].delta_position[1], -_follower_msg[index].delta_position[2]};
            const Eigen::Vector3d velocity_now = {_follower_msg[index].velocity_uav[0], _follower_msg[index].velocity_uav[1], _follower_msg[index].velocity_uav[2]};

            position_err = position_sp - position_now;
            acceleration_sp = KP * position_err - KD * velocity_now;

            RCLCPP_DEBUG(this->get_logger(), "UAV%d perr: %f %f %f", index, position_err[0], position_err[1], position_err[2]);
        };

        Eigen::Vector3d acceleration_sp, position_err[3];
        #define FILL_OUTPUT_FORCE(i) \
        { \
            cable_end_controller(i - 1, position_err[i - 1], acceleration_sp); \
            output.force_sp##i[0] = _uav_mass * acceleration_sp[0]; \
            output.force_sp##i[1] = _uav_mass * acceleration_sp[1]; \
            output.force_sp##i[2] = _uav_mass * acceleration_sp[2] - _uav_mass * 9.81; \
        }
        
        FILL_OUTPUT_FORCE(1); // leader
        FILL_OUTPUT_FORCE(2); // follower 1
        FILL_OUTPUT_FORCE(3); // follower 2

        #undef FILL_OUTPUT_FORCE
    
        // check if all uav are in position
        for (int i = 0; i < 3; i++) {
            if (position_err[i].norm() > 0.3) {
                RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, 
                        "UAV%d not in position, dist: %f", i + 1, position_err[i].norm());
                return false;
            }
        }
        return true;
    }

    void state_running(CodeGenController::OutputBus &output)
    {
        // trajectory generator
        TrajectoryGenerator::traj_out traj_out;
        _traj_gen->update(20_ms, traj_out);

        for (int i = 0; i < 3; i++) {
            _control_input.Payload_Out.pL[i] = _follower_msg[0].position_uav[i] + _follower_msg[0].delta_position[i];
            _control_input.Payload_Out.vL[i] = _follower_msg[0].velocity_uav[i] + _follower_msg[0].delta_velocity[i];
            _control_input.Payload_Out.p_1[i] = _control_input.Payload_Out.pL[i] - _follower_msg[0].delta_position[i];
            _control_input.Payload_Out.p_2[i] = _control_input.Payload_Out.pL[i] - _follower_msg[1].delta_position[i];
            _control_input.Payload_Out.p_3[i] = _control_input.Payload_Out.pL[i] - _follower_msg[2].delta_position[i];
            _control_input.Payload_Out.v_1[i] = _follower_msg[0].velocity_uav[i];
            _control_input.Payload_Out.v_2[i] = _follower_msg[1].velocity_uav[i];
            _control_input.Payload_Out.v_3[i] = _follower_msg[2].velocity_uav[i];
            _control_input.Traj_sp.pos_sp[i] = traj_out.position[i];
            _control_input.Traj_sp.vel_sp[i] = traj_out.velocity[i];
            _control_input.Traj_sp.acc_ff[i] = traj_out.acceleration[i];
        }
        _control_input.Payload_Out.timestamp = _follower_msg[0].timestamp;

        // update cable direction
        _control_input.Dir_sp.timestamp = _cable_dir_timestamp;
        for (int i = 0; i < 3; i++) {
            _control_input.Dir_sp.q_sp1[i] = _cable_dir_sp[0][i];
            _control_input.Dir_sp.q_sp2[i] = _cable_dir_sp[1][i];
            _control_input.Dir_sp.q_sp3[i] = _cable_dir_sp[2][i];
        }

        // running
        _controller.step(_control_input);

        output = _controller.getOutput();

        // test: every 5s, request force optimization
        if (1)
        {
            _disturbance_filter.update({output.state.dL[0], output.state.dL[1], output.state.dL[2]});

            static uint64_t last_request_time = _running_time;
            static Eigen::Vector3d last_trim_acc = {0, 0, -9.8};

            Eigen::Vector3d expected_trim_acc = {
                -_disturbance_filter.getState()[0],
                -_disturbance_filter.getState()[1],
                -_disturbance_filter.getState()[2] - 9.8
            };

            bool need_request = false;
            for (int i = 0; i < 3; i++) {
                if (std::abs(expected_trim_acc[i] - last_trim_acc[i]) > _margin_acc * 0.6) {
                    need_request = true;
                }
            }

            if (_running_time - last_request_time > 5_s && need_request) 
            {
                RCLCPP_INFO(this->get_logger(), "Expected: [%f %f %f], last: [%f %f %f]", 
                        expected_trim_acc[0], expected_trim_acc[1], expected_trim_acc[2], 
                        last_trim_acc[0], last_trim_acc[1], last_trim_acc[2]);

                async_request_force_opt(expected_trim_acc * _load_mass);
                last_request_time = _running_time;
                last_trim_acc = expected_trim_acc;
            }
        }

        RCLCPP_DEBUG(this->get_logger(), "\nPOS1: %f %f %f VEL1: %f %f %f\nPOS2: %f %f %f VEL2: %f %f %f\nPOS3: %f %f %f VEL3: %f %f %f\n", 
                    _control_input.Payload_Out.p_1[0], _control_input.Payload_Out.p_1[1], _control_input.Payload_Out.p_1[2], 
                    _control_input.Payload_Out.v_1[0], _control_input.Payload_Out.v_1[1], _control_input.Payload_Out.v_1[2],
                    _control_input.Payload_Out.p_2[0], _control_input.Payload_Out.p_2[1], _control_input.Payload_Out.p_2[2], 
                    _control_input.Payload_Out.v_2[0], _control_input.Payload_Out.v_2[1], _control_input.Payload_Out.v_2[2],
                    _control_input.Payload_Out.p_3[0], _control_input.Payload_Out.p_3[1], _control_input.Payload_Out.p_3[2], 
                    _control_input.Payload_Out.v_3[0], _control_input.Payload_Out.v_3[1], _control_input.Payload_Out.v_3[2]);

        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 500, "\nDisturbance: [%f %f %f]\nFilter: [%f %f %f]",
                    output.state.dL[0], output.state.dL[1], output.state.dL[2],
                    _disturbance_filter.getState()[0], _disturbance_filter.getState()[1], _disturbance_filter.getState()[2]);
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
        if (!_boot_time) {
            _boot_time = timestamp;
        }
        const uint64_t running_time = timestamp - _boot_time;
        _dt = running_time - _running_time;
        _running_time = running_time;

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
        
        RCLCPP_DEBUG(this->get_logger(), "\n UAV1: %f %f %f\n UAV2: %f %f %f\n UAV3: %f %f %f", 
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
        this->declare_parameter("traj_type", "line"); // line or circle

        // controller parameters
        this->declare_parameter("cable_len", 1.0);
        this->declare_parameter("load_mass", 1.0);
        this->declare_parameter("uav_mass", 1.5);
        this->declare_parameter("eso_enable", true);
        this->declare_parameter("eso_pl", std::vector<double>{3.0, 5.0, 10.0});
        this->declare_parameter("eso_vl", std::vector<double>{0.0, 20.0, 30.0});
        this->declare_parameter("min_tension", 0.0); // N
        this->declare_parameter("max_tension", 7.0); // N
        this->declare_parameter("kq", 1.0);
        this->declare_parameter("kw", 1.0);

        // Inmutable parameters
        _load_mass = this->get_parameter("load_mass").as_double();
        _uav_mass = this->get_parameter("uav_mass").as_double();
        _cable_len = this->get_parameter("cable_len").as_double();
        // Immutable control parameters
        CONTROL_PARAM.CABLE_LEN = _cable_len;
        CONTROL_PARAM.MASS_LOAD = _load_mass;
        CONTROL_PARAM.MASS_UAV = _uav_mass;

        // disturbance filter
        // sample freq: 50Hz, cutoff freq: 0.5Hz
        _disturbance_filter.setCutoffFreq(1_s / 20_ms, 0.1);
    }

    void parameter_update()
    {
        // Mutable control parameters
        CONTROL_PARAM.KP = 0.5;
        CONTROL_PARAM.KV = 1.0;
        CONTROL_PARAM.KQ = this->get_parameter("kq").as_double();
        CONTROL_PARAM.KW = this->get_parameter("kw").as_double();
        CONTROL_PARAM.KQI = 0.0;
        CONTROL_PARAM.TENSION_MIN = this->get_parameter("min_tension").as_double();
        CONTROL_PARAM.TENSION_MAX = this->get_parameter("max_tension").as_double();

        const bool eso_enable = this->get_parameter("eso_enable").as_bool();
        const std::vector<double> eso_pl = this->get_parameter("eso_pl").as_double_array();
        const std::vector<double> eso_vl = this->get_parameter("eso_vl").as_double_array();
        for (int i = 0; i < 3; i++) {
            CONTROL_PARAM.ESO_PL[i] = eso_enable ? eso_pl[i] : 0;
            CONTROL_PARAM.ESO_VL[i] = eso_enable ? eso_vl[i] : 0;
            CONTROL_PARAM.ESO_PI[i] = eso_enable ? eso_pl[i] : 0;
            CONTROL_PARAM.ESO_VI[i] = eso_enable ? eso_vl[i] : 0;
        }
    }

    void async_request_force_opt(const Eigen::Vector3d &center)
    {
        auto request = std::make_shared<mqsls::srv::ForceOpt::Request>();
        request->center = {center[0], center[1], center[2]};
        request->tension_max = this->get_parameter("max_tension").as_double();
        request->tension_min = this->get_parameter("min_tension").as_double();
        auto result = _force_opt_client->async_send_request(request, std::bind(&MqslsLeader::handle_force_opt_response, this, std::placeholders::_1));
    }

    void handle_force_opt_response(rclcpp::Client<mqsls::srv::ForceOpt>::SharedFuture future)
    {
        auto response = future.get();
        if (response) {
            for (int i = 0; i < 3; i++) {
                _cable_dir_sp[0][i] = response->q_1[i];
                _cable_dir_sp[1][i] = response->q_2[i];
                _cable_dir_sp[2][i] = response->q_3[i];
            }
            _margin_acc = response->radius / _load_mass;
            _cable_dir_timestamp++;
            RCLCPP_INFO(this->get_logger(), "Force optimization response:\nRadius: %f\nQ1: [%f %f %f]\nQ2: [%f %f %f]\nQ3: [%f %f %f]", 
                    response->radius, response->q_1[0], response->q_1[1], response->q_1[2], 
                    response->q_2[0], response->q_2[1], response->q_2[2], 
                    response->q_3[0], response->q_3[1], response->q_3[2]);
        } else {
            RCLCPP_ERROR(this->get_logger(), "Force optimization failed");
        }
    }
    
    // Subscriber
    rclcpp::Subscription<mqsls::msg::FollowerSend>::SharedPtr _leader_recv_sub[2];

    // Publisher
    rclcpp::Publisher<mqsls::msg::FollowerRecv>::SharedPtr _leader_send_pub[2];

    // client for force optimization
    rclcpp::Client<mqsls::srv::ForceOpt>::SharedPtr _force_opt_client;

    // detail
    uint64_t _dt = 0, _running_time = 0, _boot_time = 0;
    Eigen::Vector3d _cable_dir_sp[3] = {
        {0.25, 0.443, 0.866},
        {-0.5, 0, 0.866},
        {0.25, -0.443, 0.866}
    };
    uint64_t _cable_dir_timestamp = 0;
    double _margin_acc = 2; // m/s^2
    mqsls::msg::FollowerSend _follower_msg[3]; // 0: leader, 1: follower 2, 2: follower 3
    std::shared_ptr<InputSource> _input_source;
    std::shared_ptr<OutputActuator> _output_actuator;
    std::shared_ptr<EventHandler> _event_handler;

    // Inmutable parameters
    double _load_mass = 1.0; // kg
    double _uav_mass = 1.5; // kg
    double _cable_len = 1.0; // m

    // controller
    CodeGenController::InputBus _control_input;
    CodeGenController _controller;
    std::shared_ptr<TrajectoryGenerator> _traj_gen;
    AlphaFilter<Eigen::Vector3d> _disturbance_filter;

    // recorder
    void record()
    {
        MqslsDataFrame frame;
        frame.timestamp = _running_time;
        frame.load_position = {_control_input.Payload_Out.pL[0], _control_input.Payload_Out.pL[1], _control_input.Payload_Out.pL[2]};
        frame.load_velocity = {_control_input.Payload_Out.vL[0], _control_input.Payload_Out.vL[1], _control_input.Payload_Out.vL[2]};
        
        frame.uav_position[0] = {_control_input.Payload_Out.p_1[0], _control_input.Payload_Out.p_1[1], _control_input.Payload_Out.p_1[2]};
        frame.uav_velocity[0] = {_control_input.Payload_Out.v_1[0], _control_input.Payload_Out.v_1[1], _control_input.Payload_Out.v_1[2]};
        frame.uav_position[1] = {_control_input.Payload_Out.p_2[0], _control_input.Payload_Out.p_2[1], _control_input.Payload_Out.p_2[2]};
        frame.uav_velocity[1] = {_control_input.Payload_Out.v_2[0], _control_input.Payload_Out.v_2[1], _control_input.Payload_Out.v_2[2]};
        frame.uav_position[2] = {_control_input.Payload_Out.p_3[0], _control_input.Payload_Out.p_3[1], _control_input.Payload_Out.p_3[2]};
        frame.uav_velocity[2] = {_control_input.Payload_Out.v_3[0], _control_input.Payload_Out.v_3[1], _control_input.Payload_Out.v_3[2]};
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