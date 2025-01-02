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
    using OutputBus = mqsls::CodeGenForceOptimizer::OutputBus;

    MultiStartOptimizer()
    {
    }

    const OutputBus &optimize(const Eigen::Vector3d &center, const Eigen::Vector3d &T_min, const Eigen::Vector3d &T_max)
    {
        // Define Input
        auto input = mqsls::CodeGenForceOptimizer::InputBus();
        memcpy(input.center, center.data(), sizeof(input.center));
        memcpy(input.T_min, T_min.data(), sizeof(input.T_min));
        memcpy(input.T_max, T_max.data(), sizeof(input.T_max));

        // Optimize forces
        _best_output = _optimizer.optimize(input);

        if (_best_output.exitflag == 1 || _best_output.exitflag == 2)
        {
            std::cout << "Feasible solution found....." << std::endl;
        }
        else
        {
            std::cout << "No feasible solution found..... " << "Exitflag: " << _best_output.exitflag << std::endl;
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
    mqsls::CodeGenForceOptimizer _optimizer;
    OutputBus _best_output;
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
        auto output = _ms_opt.optimize(Eigen::Vector3d(request->center[0], request->center[1], request->center[2]), request->tension_min * Eigen::Vector3d(1, 1, 1), request->tension_max * Eigen::Vector3d(1, 1, 1));

        response->radius = output.radius;
        Eigen::Matrix3d best_result = Eigen::Map<Eigen::Matrix3d>(output.result);
        for (int i = 0; i < 3; i++)
        {
            response->q_1[i] = best_result(i, 0);
            response->q_2[i] = best_result(i, 1);
            response->q_3[i] = best_result(i, 2);
        }
    }

    MultiStartOptimizer _ms_opt;
    rclcpp::Service<mqsls::srv::ForceOpt>::SharedPtr _service;
};


int main(int argc, char **argv)
{
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<ForcePlanner>());

    rclcpp::shutdown();
    
    return 0;
}