#pragma once

extern "C"
{
    #include "force_opt/mso_forces.h"
}

namespace mqsls {

class CodeGenForceOptimizer
{
public:
    struct InputBus
    {
        double center[3];
        double T_min[3];
        double T_max[3];
    };
    
    struct OutputBus
    {
        double result[9];
        double radius;
        double exitflag;
    };
    CodeGenForceOptimizer()
    {
        mso_forces_initialize();
    }

    ~CodeGenForceOptimizer()
    {
        mso_forces_terminate();
    }

    const OutputBus &optimize(const InputBus &input)
    {
        mso_forces(input.center, input.T_min, input.T_max, _output.result, &_output.radius, &_output.exitflag);
        return _output;
    }
private:
    OutputBus _output;
};


} // namespace mqsls