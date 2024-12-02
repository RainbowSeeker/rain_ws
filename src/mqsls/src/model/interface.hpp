#pragma once

#include "controller/payload_controller.h"

namespace mqsls {

class Controller
{
public:
    using InputBus = payload_controller::ExtU_payload_controller_T;
    using OutputBus = payload_controller::ExtY_payload_controller_T;
   
    Controller()
    {
        _ctl.initialize();
    }

    ~Controller()
    {
        _ctl.terminate();
    }

    void step(InputBus &input)
    {
        _ctl.setExternalInputs(&input);
        _ctl.step();
    }

    const OutputBus &getOutput() const
    {
        return _ctl.getExternalOutputs();
    }

private:
    payload_controller _ctl;
};
    
} // namespace mqsls