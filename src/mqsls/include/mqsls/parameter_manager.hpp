#include <rclcpp/rclcpp.hpp>
#include <unordered_map>
#include <string>

namespace formation {

// Only support SingleThreadedExecutor
class ParameterManager
{
public:
    static ParameterManager & attach(rclcpp::Node *node)
    {
        ParameterManager &instance = get_instance();
        instance._node = node;
        return instance;
    }

    static int add_parameter(rclcpp::Parameter *param)
    {
        ParameterManager &instance = get_instance();
        instance._params[param->get_name()] = param;
        return 0;
    }

    int declare_parameters()
    {
        if (!_node)
            throw std::runtime_error("Node is not attached to ParameterManager");
    
        for (auto &param : _params)
        {
            _node->declare_parameter(param.first, param.second->get_parameter_value());
        }

        for (auto &param : _params)
        {
            *param.second = _node->get_parameter(param.first);
        }

        auto param_change_callback = [this](const std::vector<rclcpp::Parameter> &params) {
            auto result = rcl_interfaces::msg::SetParametersResult();
            result.successful = true;
            for (auto &param : params)
            {
                auto it = _params.find(param.get_name());
                if (it != _params.end())
                {
                    *it->second = param;
                }
            }
            return result;
        };

        _callback_handler = _node->add_on_set_parameters_callback(param_change_callback);

        return 0;
    }


private:
    ParameterManager() = default;
    ~ParameterManager() = default;

    static ParameterManager & get_instance()
    {
        static ParameterManager instance;
        return instance;
    }

    std::unordered_map<std::string, rclcpp::Parameter *> _params;
    rclcpp::Node *_node;
    rclcpp::node_interfaces::OnSetParametersCallbackHandle::SharedPtr _callback_handler;
};



} // namespace formation