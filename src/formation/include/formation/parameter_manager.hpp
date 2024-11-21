#ifndef PARAMETER_MANAGER_HPP
#define PARAMETER_MANAGER_HPP

#include <rclcpp/rclcpp.hpp>
#include <unordered_map>
#include <string>
#include <mutex>

namespace Parameter
{
RCLCPP_SMART_PTR_DEFINITIONS(rclcpp::Parameter);

// Only support SingleThreadedExecutor
class ParameterManager
{
public:
    explicit ParameterManager(rclcpp::Node *node) : _node(node)
    {
        auto param_change_callback = [this](const std::vector<rclcpp::Parameter> &params) {
            auto result = rcl_interfaces::msg::SetParametersResult();
            result.successful = true;
            for (auto &param : params)
            {
                std::lock_guard<std::mutex> lock(_mutex);
                auto it = _params.find(param.get_name());
                if (it != _params.end())
                {
                    *it->second = param;
                }
            }
            return result;
        };

        _callback_handler = _node->add_on_set_parameters_callback(param_change_callback);
    }

    ~ParameterManager() = default;

    RCUTILS_WARN_UNUSED
    Parameter::SharedPtr add_parameter(const std::string &name, const rclcpp::ParameterValue &value)
    {
        Parameter::SharedPtr param = get_parameter(name);
        if (param)
        {
            return param;
        }
    
        // create a new parameter
        param = Parameter::make_shared(name);
        // declare & get new value of parameter
        _node->declare_parameter(name, value);
        *param = _node->get_parameter(name);

        // add to map
        insert(param);
        
        return param;
    }

    template<typename ValueTypeT>
    Parameter::SharedPtr add_parameter(const std::string & name, ValueTypeT value)
    {
        return add_parameter(name, rclcpp::ParameterValue(value));
    }

    RCUTILS_WARN_UNUSED
    Parameter::SharedPtr get_parameter(const std::string &name)
    {
        std::lock_guard<std::mutex> lock(_mutex);
        auto it = _params.find(name);
        if (it != _params.end())
        {
            return it->second;
        }
        return nullptr;
    }

private:
    inline void insert(Parameter::SharedPtr param)
    {
        std::lock_guard<std::mutex> lock(_mutex);
        _params[param->get_name()] = param;
    }

    rclcpp::Node *_node;
    std::mutex _mutex;
    std::unordered_map<std::string, Parameter::SharedPtr> _params;
    rclcpp::node_interfaces::OnSetParametersCallbackHandle::SharedPtr _callback_handler;
};

} // namespace formation

#endif  // PARAMETER_MANAGER_HPP