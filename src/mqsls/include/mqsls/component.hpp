#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include "mqsls/msg/follower_send.hpp"
#include "mqsls/msg/follower_recv.hpp"

namespace mqsls
{
class InputSource
{
public:
    struct InputData
    {
        mqsls::msg::FollowerSend msg;
    };
    
    virtual void register_update_callback(std::function<void(const InputData &)> callback) = 0;
};

class OutputActuator
{
public:
    struct OutputData
    {
        mqsls::msg::FollowerRecv msg;
    };

    virtual void apply(const OutputData &) = 0;
};

class EventHandler
{
public:
    virtual void register_periodic_callback(std::function<void(uint64_t)> callback, uint64_t period_us) = 0;
};
} // namespace mqsls

#endif // !COMPONENT_HPP

