#include "3dof_mqsls_controller.hpp"

using namespace std::chrono_literals;

namespace mqsls {

#define CONTROL_PERIOD      (CONTROL_EXPORT.period * 1_ms)
#define PERIODIC_RUN(__now, __period, __func) \
        do { \
            static auto __last_time = __now; \
            if (((__now)) - __last_time > ((__period))) { \
                __func; \
                __last_time = __now; \
            } \
        } while (0)

#define DELAYED_RUN_IF(__now, __delay, __func, __cond) \
        do { \
            static auto __last_time = __now; \
            if (((__now)) - __last_time > ((__delay)) && ((__cond))) { \
                __func; \
                __last_time = __now; \
            } \
        } while (0)

template <typename T>
inline Eigen::Vector3d array_to_vector3(const T &d) {
    return Eigen::Vector3d(d[0], d[1], d[2]);
}

class MqslsController : public rclcpp::Node
{
public:
    MqslsController() : rclcpp::Node("MqslsController")
    {
        parameter_declare();
        
        rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
        auto qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);

        for (int i = 0; i < 3; i++) {
            _follower_send_sub[i] = this->create_subscription<mqsls::msg::FollowerSend>(
                "follower_send" + std::to_string(i + 1), qos,
                [this, i](const mqsls::msg::FollowerSend::SharedPtr msg) {
                    _follower_msg[i] = *msg;
                });
            
            _follower_recv_pub[i] = this->create_publisher<mqsls::msg::FollowerRecv>(
                "follower_recv" + std::to_string(i + 1), 10);
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
            if (state_running(output)) {
                // _state = STATE_STOPPING;
            }
            break;
        case STATE_STOPPING:
            if (state_stopping(output)) {
                _state = STATE_OVER;
            }
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
        auto pos_pid_controller = [&, this](int index, Eigen::Vector3d &position_err, Eigen::Vector3d &acceleration_sp) -> void
        {
            #define SINGLE_TEST 0
#if SINGLE_TEST
            TrajectoryGenerator::traj_out traj_out;
            static std::shared_ptr<TrajectoryGenerator> traj_gen {std::make_shared<CircleTrajectoryGenerator>(Eigen::Vector3d(0, 0, -0.5), 1, deg2rad(4))};
            traj_gen->update(CONTROL_PERIOD, traj_out);

            static uint64_t start_time = _running_time;
            if (_running_time - start_time < 1000_s) {
                traj_out.position = Eigen::Vector3d(0, 0, -0.5);
            }

            const Eigen::Vector3d position_sp = traj_out.position;
            const Eigen::Vector3d position_now = array_to_vector3(_follower_msg[index].position_uav);
#else
            // NED frame, payload position is the origin
            const Eigen::Vector3d q_sp = _cable_dir_sp[index];
            const double safety_len = _cable_len * 0.85; // 85% of cable length
            const Eigen::Vector3d position_sp = -q_sp * safety_len;
            const Eigen::Vector3d position_now = array_to_vector3(_follower_msg[index].position_uav) - array_to_vector3(_follower_msg[index].position_load);
#endif
            #undef SINGLE_TEST
            
            const Eigen::Vector3d velocity_now = array_to_vector3(_follower_msg[index].velocity_uav);

            position_err = position_sp - position_now;

            acceleration_sp = _pid_controller[index].computeCommand(position_err, -velocity_now, _dt);

            RCLCPP_INFO(this->get_logger(), "UAV%d pos: %f %f %f, pos_sp: %f %f %f\nacc_sp: %f %f %f", 
                        index, position_now[0], position_now[1], position_now[2], 
                        position_sp[0], position_sp[1], position_sp[2], 
                        acceleration_sp[0], acceleration_sp[1], acceleration_sp[2]);
            RCLCPP_DEBUG(this->get_logger(), "UAV%d perr: %f %f %f", index, position_err[0], position_err[1], position_err[2]);
        };

        Eigen::Vector3d acceleration_sp, position_err[3];
        #define FILL_OUTPUT_FORCE(i) \
        { \
            pos_pid_controller(i - 1, position_err[i - 1], acceleration_sp); \
            output.force_sp##i[0] = _uav_mass * acceleration_sp[0]; \
            output.force_sp##i[1] = _uav_mass * acceleration_sp[1]; \
            output.force_sp##i[2] = _uav_mass * acceleration_sp[2] - _uav_mass * 9.81; \
        }
        
        FILL_OUTPUT_FORCE(1); // leader
        FILL_OUTPUT_FORCE(2); // follower 1
        FILL_OUTPUT_FORCE(3); // follower 2

        #undef FILL_OUTPUT_FORCE
    
        // check if all uav are in position
        bool all_in_position = true;
        for (int i = 0; i < 3; i++) {
            if (position_err[i].norm() > 0.3) {
                RCLCPP_INFO(this->get_logger(), "UAV%d not in position, dist: %f", i + 1, position_err[i].norm());
                all_in_position = false;
            }
            else {
                RCLCPP_INFO(this->get_logger(), "UAV%d READY, dist: %f", i + 1, position_err[i].norm());
            }
        }

        return all_in_position;
    }

    bool state_running(CodeGenController::OutputBus &output)
    {
        // trajectory generator
        TrajectoryGenerator::traj_out traj_out;
        _traj_gen->update(CONTROL_PERIOD, traj_out);

        // controller step
        output = controller_step(traj_out);

        // ACCS : every 5s, request force optimization
        auto accs_update = [this, &output]() -> void
        {
            static AlphaFilter<Eigen::Vector3d> trim_acc_filter(0.5, {0, 0, -9.8});

            Eigen::Vector3d expected_trim_acc = {
                -output.state.dL[0],
                -output.state.dL[1],
                -output.state.dL[2] - 9.8
            };
            trim_acc_filter.update(expected_trim_acc);

            async_request_force_opt(trim_acc_filter.getState() * _load_mass);

            RCLCPP_INFO(this->get_logger(), "Request ForceOpt: [%f %f %f]", trim_acc_filter.getState()[0], trim_acc_filter.getState()[1], trim_acc_filter.getState()[2]);
        };

        const bool accs_enable = this->get_parameter("accs_enable").as_bool();
        if (accs_enable)
        {
#if 1
            #define EPSILON 1e-3
            static int sp_count = 0;
            if (_traj_type == "line" && std::abs(traj_out.position.x() - 0) < EPSILON) {
                sp_count++;
                if (sp_count == 1) {
                    accs_update();
                }
            }
#else
            bool need_accs_update = output.state.margin < _opt_margin * 0.2; // 20% threshold
            DELAYED_RUN_IF(_running_time, 
                            5_s, 
                            accs_update(),
                            need_accs_update
            );
#endif
        }
        RCLCPP_DEBUG(this->get_logger(), "\nPOS1: %f %f %f VEL1: %f %f %f\nPOS2: %f %f %f VEL2: %f %f %f\nPOS3: %f %f %f VEL3: %f %f %f\n", 
                    _control_input.Payload_Out.p_1[0], _control_input.Payload_Out.p_1[1], _control_input.Payload_Out.p_1[2], 
                    _control_input.Payload_Out.v_1[0], _control_input.Payload_Out.v_1[1], _control_input.Payload_Out.v_1[2],
                    _control_input.Payload_Out.p_2[0], _control_input.Payload_Out.p_2[1], _control_input.Payload_Out.p_2[2], 
                    _control_input.Payload_Out.v_2[0], _control_input.Payload_Out.v_2[1], _control_input.Payload_Out.v_2[2],
                    _control_input.Payload_Out.p_3[0], _control_input.Payload_Out.p_3[1], _control_input.Payload_Out.p_3[2], 
                    _control_input.Payload_Out.v_3[0], _control_input.Payload_Out.v_3[1], _control_input.Payload_Out.v_3[2]);

        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 500, "Disturbance: [%.2f %.2f %.2f], Margin: %.2f%%",output.state.dL[0], output.state.dL[1], output.state.dL[2], output.state.margin * 100 / _opt_margin);
        RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, "running time: %f", (double)_running_time / 1_s);
        
        if (_running_time < _lasting_time) {
            return false;
        }
        return true;
    }

    bool state_stopping(CodeGenController::OutputBus &output)
    {
        // stopping

        // move to landing position
        static const TrajectoryGenerator::traj_out traj_out = {
            {_follower_msg[0].position_load[0], _follower_msg[0].position_load[1], 0},  // position
            {0, 0, 0},  // velocity
            {0, 0, 0}   // acceleration
        };

        // controller step
        output = controller_step(traj_out);

        // check if load is in position
        Eigen::Vector3d position_err = array_to_vector3(_follower_msg[0].position_load) - traj_out.position;
        if (position_err.norm() > 0.3) {
            RCLCPP_INFO_THROTTLE(this->get_logger(), *get_clock(), 1000, 
                    "Load not in position, dist: %f", position_err.norm());
            return false;
        }
        return true;
    }

    void state_over(CodeGenController::OutputBus &output)
    {
        // over
        // uav landing with 0.5 m/s velocity
        const double KP = 1;
        const double KD = 1.2;
        const double landing_speed = 0.5;

        #define FILL_OUTPUT_FORCE(i) \
        { \
            static const Eigen::Vector3d position_sp = array_to_vector3(_follower_msg[i - 1].position_uav); \
            const Eigen::Vector3d position_now = array_to_vector3(_follower_msg[i - 1].position_uav); \
            const Eigen::Vector3d velocity_now = array_to_vector3(_follower_msg[i - 1].velocity_uav); \
            const Eigen::Vector3d position_err = (position_sp - position_now).cwiseProduct(Eigen::Vector3d(1, 1, 0)) + landing_speed * Eigen::Vector3d::UnitZ(); \
            Eigen::Vector3d acceleration_sp = KP * position_err - KD * velocity_now; \
            output.force_sp##i[0] = _uav_mass * acceleration_sp[0]; \
            output.force_sp##i[1] = _uav_mass * acceleration_sp[1]; \
            output.force_sp##i[2] = _uav_mass * acceleration_sp[2] - _uav_mass * 9.81; \
        }
        
        FILL_OUTPUT_FORCE(1); // leader
        FILL_OUTPUT_FORCE(2); // follower 1
        FILL_OUTPUT_FORCE(3); // follower 2

        #undef FILL_OUTPUT_FORCE
    }

    void input_update()
    {
        for (int i = 0; i < 3; i++) {
            _control_input.Payload_Out.pL[i] = _follower_msg[0].position_load[i];
            _control_input.Payload_Out.vL[i] = _follower_msg[0].velocity_load[i];
            _control_input.Payload_Out.p_1[i] = _control_input.Payload_Out.pL[i] + _follower_msg[0].position_uav[i] - _follower_msg[0].position_load[i];
            _control_input.Payload_Out.p_2[i] = _control_input.Payload_Out.pL[i] + _follower_msg[1].position_uav[i] - _follower_msg[1].position_load[i];
            _control_input.Payload_Out.p_3[i] = _control_input.Payload_Out.pL[i] + _follower_msg[2].position_uav[i] - _follower_msg[2].position_load[i];
            _control_input.Payload_Out.v_1[i] = _follower_msg[0].velocity_uav[i];
            _control_input.Payload_Out.v_2[i] = _follower_msg[1].velocity_uav[i];
            _control_input.Payload_Out.v_3[i] = _follower_msg[2].velocity_uav[i];
        }
        _control_input.Payload_Out.timestamp = _follower_msg[0].timestamp;

        // update cable direction
        _control_input.Dir_sp.timestamp = _cable_dir_timestamp;
        for (int i = 0; i < 3; i++) {
            _control_input.Dir_sp.q_sp1[i] = _cable_dir_sp[0][i];
            _control_input.Dir_sp.q_sp2[i] = _cable_dir_sp[1][i];
            _control_input.Dir_sp.q_sp3[i] = _cable_dir_sp[2][i];
        }
    }

    const CodeGenController::OutputBus &controller_step(const TrajectoryGenerator::traj_out &traj_out)
    {
        for (int i = 0; i < 3; i++) {
            _control_input.Traj_sp.pos_sp[i] = traj_out.position[i];
            _control_input.Traj_sp.vel_sp[i] = traj_out.velocity[i];
            _control_input.Traj_sp.acc_ff[i] = traj_out.acceleration[i];
        }

        // running
        _controller.step(_control_input);

        return _controller.getOutput();
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

        input_update();

        CodeGenController::OutputBus output;

        state_machine_flow(output);
        
        #define FILL_OUTPUT_ACTUATOR(i) \
        { \
            auto msg = mqsls::msg::FollowerRecv(); \
            msg.timestamp = _running_time; \
            msg.force[0] = output.force_sp##i[0]; \
            msg.force[1] = output.force_sp##i[1]; \
            msg.force[2] = output.force_sp##i[2]; \
            msg.thrust_coff = _hover_thrust / _uav_mass / 9.81; \
            _follower_recv_pub[i - 1]->publish(msg); \
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
        // parameters check
        if (_test_phase != "sil" && _test_phase != "hil" && _test_phase != "actual") {
            throw std::runtime_error("Invalid test phase: " + _test_phase);
        }

        // controller parameters
        this->declare_parameter("eso_enable", true);
        this->declare_parameter("pwas_enable", true);
        this->declare_parameter("accs_enable", true);
        this->declare_parameter("eso_pl", std::vector<double>{3.0, 5.0, 10.0});
        this->declare_parameter("eso_vl", std::vector<double>{0.0, 20.0, 30.0});
        this->declare_parameter("min_tension", 0.0); // N
        this->declare_parameter("max_tension", 7.0); // N
        this->declare_parameter("kp", 0.5);
        this->declare_parameter("kv", 1.0);
        this->declare_parameter("kq", 1.0);
        this->declare_parameter("kw", 1.0);

        this->declare_parameter("traj_type", "point");

        // Immutable control parameters
        CONTROL_PARAM.CABLE_LEN = _cable_len;
        CONTROL_PARAM.MASS_LOAD = _load_mass;
        CONTROL_PARAM.MASS_UAV = _uav_mass;

        const double MC_POS_KP = 3.4;
        const double MC_POS_KI = 0.2;
        const double MC_POS_KD = 2.0;
        const double MC_POS_KI_MAX = 3.0;
        const double MC_POS_KI_MIN = -3.0;
        for (int i = 0; i < 3; i++) {
            _pid_controller[i].initPid(
                    MC_POS_KP * Eigen::Vector3d::Ones(), 
                    MC_POS_KI * Eigen::Vector3d::Ones(),
                    MC_POS_KD * Eigen::Vector3d::Ones(),
                    MC_POS_KI_MAX * Eigen::Vector3d::Ones(),
                    MC_POS_KI_MIN * Eigen::Vector3d::Ones(),
                    true);
        }
    }

    void parameter_update()
    {
        // Mutable control parameters
        CONTROL_PARAM.KP = this->get_parameter("kp").as_double();
        CONTROL_PARAM.KV = this->get_parameter("kv").as_double();
        CONTROL_PARAM.KQ = this->get_parameter("kq").as_double();
        CONTROL_PARAM.KW = this->get_parameter("kw").as_double();

        const bool pwas_enable = this->get_parameter("pwas_enable").as_bool();
        CONTROL_PARAM.TENSION_MIN = pwas_enable ? this->get_parameter("min_tension").as_double() : -1000;
        CONTROL_PARAM.TENSION_MAX = pwas_enable ? this->get_parameter("max_tension").as_double() : 1000;

        const bool eso_enable = this->get_parameter("eso_enable").as_bool();
        const std::vector<double> eso_pl = this->get_parameter("eso_pl").as_double_array();
        const std::vector<double> eso_vl = this->get_parameter("eso_vl").as_double_array();
        for (int i = 0; i < 3; i++) {
            CONTROL_PARAM.ESO_PL[i] = eso_enable ? eso_pl[i] : 0;
            CONTROL_PARAM.ESO_VL[i] = eso_enable ? eso_vl[i] : 0;
            CONTROL_PARAM.ESO_PI[i] = eso_enable ? eso_pl[i] : 0;
            CONTROL_PARAM.ESO_VI[i] = eso_enable ? eso_vl[i] : 0;
        }
        CONTROL_PARAM.KPI = eso_enable ? 0 : 0.2;
        CONTROL_PARAM.KQI = eso_enable ? 0 : 0.2;

        const std::string new_traj_type = this->get_parameter("traj_type").as_string();
        if (new_traj_type != _traj_type) {
            auto traj_gen = make_trajectory_generator(new_traj_type, _test_phase);
            if (traj_gen) {
                _traj_gen = traj_gen;
                _traj_type = new_traj_type;
                RCLCPP_INFO(this->get_logger(), "Trajectory type changed to: %s", new_traj_type.c_str());
            }
            else {
                RCLCPP_ERROR(this->get_logger(), "Invalid trajectory type: %s", new_traj_type.c_str());
            }
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
            _opt_margin = response->radius;
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
    rclcpp::Subscription<mqsls::msg::FollowerSend>::SharedPtr _follower_send_sub[3];// 0: leader, 1: follower 2, 2: follower 3

    // Publisher
    rclcpp::Publisher<mqsls::msg::FollowerRecv>::SharedPtr _follower_recv_pub[3]; // 0: leader, 1: follower 2, 2: follower 3

    // client for force optimization
    rclcpp::Client<mqsls::srv::ForceOpt>::SharedPtr _force_opt_client;

    // Timer
    rclcpp::TimerBase::SharedPtr _timer;

    // Inmutable parameters
    const std::string _test_phase = this->declare_parameter("test_phase", "sil"); // sil, hil, actual
    const double _load_mass = this->declare_parameter("load_mass", 1.0); // kg
    const double _uav_mass = this->declare_parameter("uav_mass", 1.5); // kg
    const double _cable_len = this->declare_parameter("cable_len", 1.0); // m
    const double _hover_thrust = this->declare_parameter("hover_thrust", 0.5); // percent %
    const uint64_t _lasting_time = this->declare_parameter("lasting_time", 60) * 1_s; // s

    // detail
    uint64_t _dt = 0, _running_time = 0, _boot_time = 0;
    Eigen::Vector3d _cable_dir_sp[3] = {
        {0.25, 0.443, 0.866},
        {-0.5, 0, 0.866},
        {0.25, -0.443, 0.866}
    };
    uint64_t _cable_dir_timestamp = 0;
    double _opt_margin = 2; // N
    mqsls::msg::FollowerSend _follower_msg[3]; // 0: leader, 1: follower 2, 2: follower 3

    // controller
    CodeGenController::InputBus _control_input;
    CodeGenController _controller;
    mqsls::Pid<Eigen::Vector3d> _pid_controller[3]; // 0: leader, 1: follower 2, 2: follower 3
    std::string _traj_type {"invalid"};
    std::shared_ptr<TrajectoryGenerator> _traj_gen;

    // recorder
    void record()
    {
        MqslsDataFrame frame;
        frame.timestamp = _running_time;
        frame.load_position = array_to_vector3(_control_input.Payload_Out.pL);
        frame.load_velocity = array_to_vector3(_control_input.Payload_Out.vL);
        
        frame.uav_position[0] = array_to_vector3(_control_input.Payload_Out.p_1);
        frame.uav_velocity[0] = array_to_vector3(_control_input.Payload_Out.v_1);
        frame.uav_position[1] = array_to_vector3(_control_input.Payload_Out.p_2);
        frame.uav_velocity[1] = array_to_vector3(_control_input.Payload_Out.v_2);
        frame.uav_position[2] = array_to_vector3(_control_input.Payload_Out.p_3);
        frame.uav_velocity[2] = array_to_vector3(_control_input.Payload_Out.v_3);
        
        frame.dL = array_to_vector3(_controller.getOutput().state.dL);
        frame.tau = array_to_vector3(_controller.getOutput().state.tau);
        frame.margin = _controller.getOutput().state.margin;

        frame.load_position_sp = array_to_vector3(_control_input.Traj_sp.pos_sp);
        frame.load_velocity_sp = array_to_vector3(_control_input.Traj_sp.vel_sp);
        frame.cable_dir_sp[0] = array_to_vector3(_control_input.Dir_sp.q_sp1);
        frame.cable_dir_sp[1] = array_to_vector3(_control_input.Dir_sp.q_sp2);
        frame.cable_dir_sp[2] = array_to_vector3(_control_input.Dir_sp.q_sp3);
        
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