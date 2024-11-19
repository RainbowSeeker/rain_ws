#pragma once

#include <string>
#include <iomanip>
#include <cmath>
#include <cstdint>

#ifndef MATH_PI
#define MATH_PI		3.141592653589793238462643383280
#endif

namespace math {

template <typename T>
T constrain(T value, T min_value, T max_value) {
    if (value < min_value) {
        return min_value;
    } else if (value > max_value) {
        return max_value;
    } else {
        return value;
    }
}


template<typename _Tp>
constexpr bool isInRange(_Tp val, _Tp min_val, _Tp max_val)
{
	return (min_val <= val) && (val <= max_val);
}

template<typename T>
constexpr T radians(T degrees)
{
	return degrees * (static_cast<T>(MATH_PI) / static_cast<T>(180));
}

template<typename T>
constexpr T degrees(T radians)
{
	return radians * (static_cast<T>(180) / static_cast<T>(MATH_PI));
}

}


namespace utils {

std::string realtime_str()
{
    // format: MM-DD HH_MM_SS
    auto now = std::time(nullptr);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&now), "%m-%d %H_%M_%S");
    return ss.str();
}

}

// literals
constexpr double operator""_deg(long double deg)
{
    return math::radians(deg);
}

constexpr uint64_t operator""_us(long double us)
{
    return us;
}

constexpr uint64_t operator""_us(unsigned long long us)
{
    return us;
}

constexpr uint64_t operator""_ms(long double ms)
{
    return ms * 1e3;
}

constexpr uint64_t operator""_ms(unsigned long long ms)
{
    return ms * 1e3;
}

constexpr uint64_t operator""_s(long double s)
{
    return s * 1e6;
}

constexpr uint64_t operator""_s(unsigned long long s)
{
    return s * 1e6;
}