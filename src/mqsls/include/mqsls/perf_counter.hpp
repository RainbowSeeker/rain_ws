#pragma once

#include <string>

enum class PerfCounterType
{
    INTERVAL,
    COUNTER
};

class PerfCounter
{
public:
    PerfCounter(PerfCounterType type, const std::string &name)
        : _type(type), _name(name)
    {
    }

    uint64_t count() const
    {
        return _count;
    }

    void begin()
    {
        _start = std::chrono::high_resolution_clock::now();
        _is_running = true;
    }

    void end()
    {
        if (!_is_running)
        {
            return;
        }
        
        _end = std::chrono::high_resolution_clock::now();
        _duration = std::chrono::duration_cast<std::chrono::microseconds>(_end - _start);
        _count++;

        // update statistics
        if ((uint64_t)_duration.count() < _min)   _min = _duration.count();
        if ((uint64_t)_duration.count() > _max)   _max = _duration.count();

        _mean = (_mean * _count + _duration.count()) / (_count + 1);
        _variance = (_variance * _count + (_duration.count() - _mean) * (_duration.count() - _mean)) / (_count + 1);
        _stddev = sqrt(_variance);
    }

    void tick()
    {
        end();
        begin();
    }

    void print()
    {
        switch (_type)
        {
        case PerfCounterType::INTERVAL:
            std::cout << "Interval [" << _name << "]: " << "count: " << _count << ", " << "mean: " << _mean << "us, " << "stddev: " << _stddev << "us, " << "min: " << _min << "us, " << "max: " << _max << "us" << std::endl;
            break;
        case PerfCounterType::COUNTER:
            std::cout << "Counter [" << _name << "]: " << _count << std::endl;
            break;
        }
    }

private:
    const PerfCounterType _type;
    const std::string _name;

    bool _is_running = false;
    std::chrono::time_point<std::chrono::high_resolution_clock> _start, _end;
    std::chrono::microseconds _duration;
    uint64_t _count = 0;

    // statistics
    uint64_t _min = std::numeric_limits<uint64_t>::max(), _max = std::numeric_limits<uint64_t>::min();
    double _mean = 0;
    double _variance = 0;
    double _stddev = 0;
};