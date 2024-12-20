#include <chrono>
#include <vector>
#include <iostream>
#include <cstring>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <rclcpp/rclcpp.hpp>
#include <mqsls/srv/force_opt.hpp>

#include "../model/interface.hpp"

#define deg2rad(x) ((x) * M_PI / 180)
#define rad2deg(x) ((x) * 180 / M_PI)
#define TIMEIT(code) \
        { \
            auto __start = std::chrono::high_resolution_clock::now(); \
            code; \
            auto __end = std::chrono::high_resolution_clock::now(); \
            std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count() << "us" << std::endl; \
        }

class MultiStartOptimizer
{
public:
    struct Configuration
    {
        std::vector<Eigen::Vector3d> theta;
        std::vector<Eigen::Vector3d> psi;
    };
    static Configuration default_config;

    using OutputBus = mqsls::CodeGenForceOptimizer::OutputBus;

    MultiStartOptimizer(const Configuration &config = default_config, bool early_exit = true) : _early_exit(early_exit)
    {
        _all_initial_guesses.reserve(config.theta.size() * config.psi.size());
        for (const auto &theta : config.theta)
        {
            for (const auto &psi : config.psi)
            {
                Eigen::Matrix3d Q;
                for (int i = 0; i < 3; i++)
                {
                    Q.col(i) = Eigen::Vector3d(-cos(psi[i]) * cos(theta[i]), -sin(psi[i]) * cos(theta[i]), sin(theta[i]));
                }
                _all_initial_guesses.push_back(Q);
            }
        }
    }

    const OutputBus &optimize(const Eigen::Vector3d &center, const double T_max, const Eigen::Vector3d &psi = {NAN, NAN, NAN})
    {
        _best_output = _suboptimal_output = OutputBus();

        const Eigen::Vector3d T_min_vec = {0.0, 0.0, 0.0};
        const Eigen::Vector3d T_max_vec = {T_max, T_max, T_max};

        // Define Input
        auto input = mqsls::CodeGenForceOptimizer::InputBus();
        memcpy(input.center, center.data(), sizeof(input.center));
        memcpy(input.T_min, T_min_vec.data(), sizeof(input.T_min));
        memcpy(input.T_max, T_max_vec.data(), sizeof(input.T_max));
        memcpy(input.psi, psi.data(), sizeof(input.psi));

        // Optimize forces for each initial guess
        int iter = 1;
        for (const auto &initial_guess : _all_initial_guesses)
        {
            memcpy(input.initial_guess, initial_guess.data(), sizeof(input.initial_guess));
            auto output = _optimizer.optimize(input);

            if (output.exitflag == 1 || output.exitflag == 2)
            {
                if (output.radius > _best_output.radius)
                {
                    _best_output = output;
                }
                else if (_early_exit)
                {
                    std::cout << "Early exit at iteration [" << iter << "]....." << std::endl;
                    break;
                }
            }
            else
            {
                if (output.radius > _suboptimal_output.radius)
                {
                    _suboptimal_output = output;
                }
            }
            iter++;
        }
        if (_best_output.radius < _suboptimal_output.radius)
        {
            _best_output = _suboptimal_output;
            std::cout << "No feasible solution found, using suboptimal solution....." << std::endl;
        }
        else
        {
            std::cout << "Feasible solution found....." << std::endl;
        }

        Eigen::Matrix3d best_result = Eigen::Map<Eigen::Matrix3d>(_best_output.result);
        for (int i = 0; i < 3; i++)
        {
            best_result.col(i).normalize();
        }
        memcpy(_best_output.result, best_result.data(), sizeof(_best_output.result));

        // Print output
        Eigen::Vector3d best_theta, best_psi;
        for (int i = 0; i < 3; i++)
        {
            best_theta[i] = rad2deg(atan2(best_result(2, i), sqrt(best_result(0, i) * best_result(0, i) + best_result(1, i) * best_result(1, i))));
            best_psi[i] = rad2deg(atan2(-best_result(1, i), -best_result(0, i)));
        }

        // set precision
        auto old = std::cout.precision(4);
        std::cout << "Exitflag: \t" << _best_output.exitflag << std::endl;
        std::cout << "Ball center: \t[" << center.transpose() << "]" << std::endl;
        std::cout << "Ball Radius: \t" << _best_output.radius << std::endl;
        std::cout << "Result: \t-------------------" << std::endl;
        std::cout << "  theta: \t[" << best_theta.transpose() << "]" << std::endl;
        std::cout << "  psi: \t\t[" << best_psi.transpose() << "]" << std::endl;
        // unset precision
        std::cout.precision(old);

        return _best_output;
    }

private:
    // Parameters
    const bool _early_exit = false;

    std::vector<Eigen::Matrix3d> _all_initial_guesses;
    mqsls::CodeGenForceOptimizer _optimizer;
    OutputBus _best_output, _suboptimal_output;
};

MultiStartOptimizer::Configuration MultiStartOptimizer::default_config = {
    {
        {deg2rad(10), deg2rad(10), deg2rad(10)},
        {deg2rad(20), deg2rad(20), deg2rad(20)},
        {deg2rad(30), deg2rad(30), deg2rad(30)},
        {deg2rad(40), deg2rad(40), deg2rad(40)},
        {deg2rad(50), deg2rad(50), deg2rad(50)},
        {deg2rad(60), deg2rad(60), deg2rad(60)},
        {deg2rad(70), deg2rad(70), deg2rad(70)},
        {deg2rad(80), deg2rad(80), deg2rad(80)},
    },
    {
        {deg2rad(-120), deg2rad(0), deg2rad(120)},
    }
};

class ForcePlanner : public rclcpp::Node
{
public:
    ForcePlanner() : Node("force_planner")
    {
        // Service
        _service = this->create_service<mqsls::srv::ForceOpt>("force_opt", std::bind(&ForcePlanner::handle_force_opt, this, std::placeholders::_1, std::placeholders::_2));

        RCLCPP_INFO(this->get_logger(), "Force planner ready.....");
    }

private:
    void handle_force_opt(const std::shared_ptr<mqsls::srv::ForceOpt::Request> request, std::shared_ptr<mqsls::srv::ForceOpt::Response> response)
    {
        auto output = _optimizer.optimize(Eigen::Vector3d(request->center[0], request->center[1], request->center[2]), request->tension_max, Eigen::Vector3d(request->psi[0], request->psi[1], request->psi[2]));

        response->radius = output.radius;
        Eigen::Matrix3d best_result = Eigen::Map<Eigen::Matrix3d>(output.result);
        for (int i = 0; i < 3; i++)
        {
            response->q_1[i] = best_result(i, 0);
            response->q_2[i] = best_result(i, 1);
            response->q_3[i] = best_result(i, 2);
        }
    }

    MultiStartOptimizer _optimizer;
    rclcpp::Service<mqsls::srv::ForceOpt>::SharedPtr _service;
};


int main(int argc, char **argv)
{
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<ForcePlanner>());

    rclcpp::shutdown();
    
    return 0;
}