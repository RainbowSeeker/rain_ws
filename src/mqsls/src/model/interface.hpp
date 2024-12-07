#pragma once

#include "controller/control_3dof.h"

namespace mqsls {

class CodeGenController
{
#if 1
public:
    using InputBus = control_3dof::ExtU_control_3dof_T;
    using OutputBus = control_3dof::ExtY_control_3dof_T;
   
    CodeGenController()
    {
        _ctl.initialize();
    }

    ~CodeGenController()
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
    control_3dof _ctl;
#else
public:
    using InputBus = ExtU_control_3dof_T;
    using OutputBus = ExtY_control_3dof_T;
   
    CodeGenController()
    {
        control_3dof_initialize();
    }

    ~CodeGenController()
    {
        control_3dof_terminate();
    }

    void step(InputBus &input)
    {
        control_3dof_U = input;
        control_3dof_step();
    }

    const OutputBus &getOutput() const
    {
        return control_3dof_Y;
    }
#endif
};
    
} // namespace mqsls