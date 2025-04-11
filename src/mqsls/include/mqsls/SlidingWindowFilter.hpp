#pragma once

#include <vector>
#include <stdexcept>

template <typename T>
class SlidingWindowFilter
{
public:
    using vector = std::vector<T>;

    explicit SlidingWindowFilter(int window_size) : SlidingWindowFilter(window_size, T()) {}

    explicit SlidingWindowFilter(int window_size, const T &initial_state) : _window_size(window_size), _filter_state(initial_state) 
    {
        if (window_size <= 0) {
            throw std::invalid_argument("Window size must be greater than zero.");
        }

        _samples.reserve(window_size);
    }

    ~SlidingWindowFilter() = default;

    const T &update(const T &sample)
    {
        if ((int)_samples.size() < _window_size) {
            _samples.push_back(sample);
        } else {
            _sum = _sum - _samples[_index];
            _samples[_index] = sample;
        }
        _sum = _sum + sample;

        _filter_state = _sum / _samples.size();
        _index = (_index + 1) % _window_size;

        return _filter_state;
    }

    const T &getState() const
    {
        return _filter_state;
    }
private:
    int _window_size, _index = 0;
    vector _samples;
    T _filter_state {}, _sum {};
};


