#include "3dof_mqsls_controller.hpp"

using namespace std::chrono_literals;

namespace mqsls {

#define CONTROL_PERIOD      (CONTROL_EXPORT.period * 1_ms)

class MqslsController : public rclcpp::Node
{
public:
    MqslsController() : rclcpp::Node("MqslsController")
    {
        parameter_declare();

        const std::string traj_type = this->get_parameter("traj_type").as_string();
        
        _traj_gen = make_trajectory_generator(traj_type);

        rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
        auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);

        for (int i = 0; i < 3; i++) {
            _follower_send_sub[i] = this->create_subscription<mqsls::msg::FollowerSend>(
                "follower_send" + std::to_string(i + 1), qos,
                [this, i](const mqsls::msg::FollowerSend::SharedPtr msg) {
                    _follower_msg[i] = *msg;
                });
            
            _output_actuator[i] = std::make_shared<PX4OutputActuator>(this, i + 1);
        }

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
        _timer = this->create_wall_timer(std::chrono::microseconds(CONTROL_PERIOD), std::bind(&MqslsController::run, this));
    }
    ~MqslsController() = default;
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
            const Eigen::Vector3d position_now = {
                _follower_msg[index].position_uav[0] - _follower_msg[index].position_load[0],
                _follower_msg[index].position_uav[1] - _follower_msg[index].position_load[1],
                _follower_msg[index].position_uav[2] - _follower_msg[index].position_load[2],
            };
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
        _traj_gen->update(CONTROL_PERIOD, traj_out);

        for (int i = 0; i < 3; i++) {
            _control_input.Payload_Out.pL[i] = _follower_msg[0].position_load[i];
            _control_input.Payload_Out.vL[i] = _follower_msg[0].velocity_load[i];
            _control_input.Payload_Out.p_1[i] = _control_input.Payload_Out.pL[i] + _follower_msg[0].position_uav[i] - _follower_msg[0].position_load[i];
            _control_input.Payload_Out.p_2[i] = _control_input.Payload_Out.pL[i] + _follower_msg[1].position_uav[i] - _follower_msg[1].position_load[i];
            _control_input.Payload_Out.p_3[i] = _control_input.Payload_Out.pL[i] + _follower_msg[2].position_uav[i] - _follower_msg[2].position_load[i];
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

            bool need_request = output.state.margin < _margin * 0.2; // 20% threshold

            if (need_request && _running_time - last_request_time > 5_s) 
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

        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 500, "Disturbance: [%f %f %f] \tFilter: [%f %f %f]",
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

    void run()
    {
        const uint64_t timestamp = this->get_clock()->now().nanoseconds() / 1e3;
        if (!_boot_time) {
            _boot_time = timestamp;
        }
        const uint64_t running_time = timestamp - _boot_time;
        _dt = running_time - _running_time;
        _running_time = running_time;

        parameter_update();

        CodeGenController::OutputBus output;

        state_machine_flow(output);
        
        #define FILL_OUTPUT_ACTUATOR(i) \
        { \
            auto msg = mqsls::msg::FollowerRecv(); \
            msg.timestamp = _running_time; \
            msg.force[0] = output.force_sp##i[0]; \
            msg.force[1] = output.force_sp##i[1]; \
            msg.force[2] = output.force_sp##i[2]; \
            _output_actuator[i - 1]->apply(msg); \
        }

        FILL_OUTPUT_ACTUATOR(1);
        FILL_OUTPUT_ACTUATOR(2);
        FILL_OUTPUT_ACTUATOR(3);

        #undef FILL_OUTPUT_ACTUATOR

        record();
        
        RCLCPP_DEBUG(this->get_logger(), "\n UAV1: %f %f %f\n UAV2: %f %f %f\n UAV3: %f %f %f", 
                    output.force_sp1[0], output.force_sp1[1], output.force_sp1[2],
                    output.force_sp2[0], output.force_sp2[1], output.force_sp2[2],
                    output.force_sp3[0], output.force_sp3[1], output.force_sp3[2]);
    }

    // parameters
    void parameter_declare()
    {
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
        _disturbance_filter.setCutoffFreq(1_s / CONTROL_PERIOD, 0.1);
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
        auto result = _force_opt_client->async_send_request(request, std::bind(&MqslsController::handle_force_opt_response, this, std::placeholders::_1));
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
            _margin = response->radius;
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
    rclcpp::Subscription<mqsls::msg::FollowerSend>::SharedPtr _follower_send_sub[3];

    // client for force optimization
    rclcpp::Client<mqsls::srv::ForceOpt>::SharedPtr _force_opt_client;

    // Timer
    rclcpp::TimerBase::SharedPtr _timer;

    // detail
    uint64_t _dt = 0, _running_time = 0, _boot_time = 0;
    Eigen::Vector3d _cable_dir_sp[3] = {
        {0.25, 0.443, 0.866},
        {-0.5, 0, 0.866},
        {0.25, -0.443, 0.866}
    };
    uint64_t _cable_dir_timestamp = 0;
    double _margin = 2; // N
    mqsls::msg::FollowerSend _follower_msg[3]; // 0: leader, 1: follower 2, 2: follower 3
    std::shared_ptr<PX4OutputActuator> _output_actuator[3]; // 0: leader, 1: follower 2, 2: follower 3

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
} // namespace mqsls


int main(int argc, const char** argv) 
{
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<mqsls::MqslsController>());

    rclcpp::shutdown();
    
    return 0;
}