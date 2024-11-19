#pragma once

#include <cmath>

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

