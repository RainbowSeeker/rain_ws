#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <iomanip>
#include <cmath>
#include <cstdint>
#include <Eigen/Eigen>

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

std::string nowstr()
{
    // format: MM-DD HH_MM_SS
    auto now = std::time(nullptr);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&now), "%m-%d %H_%M_%S");
    return ss.str();
}

// Quaternion
namespace quaternion
{

Eigen::Quaterniond quaternion_from_euler(const Eigen::Vector3d &euler)
{
	// YPR is ZYX axes
	return Eigen::Quaterniond(Eigen::AngleAxisd(euler.z(), Eigen::Vector3d::UnitZ()) *
				  Eigen::AngleAxisd(euler.y(), Eigen::Vector3d::UnitY()) *
				  Eigen::AngleAxisd(euler.x(), Eigen::Vector3d::UnitX()));
}

Eigen::Quaterniond quaternion_from_euler(const double roll, const double pitch, const double yaw)
{
	return quaternion_from_euler(Eigen::Vector3d(roll, pitch, yaw));
}

Eigen::Vector3d quaternion_to_euler(const Eigen::Quaterniond &q)
{
	// YPR is ZYX axes
	return q.toRotationMatrix().eulerAngles(2, 1, 0).reverse();
}

void quaternion_to_euler(const Eigen::Quaterniond &q, double &roll, double &pitch, double &yaw)
{
	const auto euler = quaternion_to_euler(q);
	roll = euler.x();
	pitch = euler.y();
	yaw = euler.z();
}

void eigen_quat_to_array(const Eigen::Quaterniond &q, std::array<float, 4> &qarray)
{
	qarray[0] = q.w();
	qarray[1] = q.x();
	qarray[2] = q.y();
	qarray[3] = q.z();
}

Eigen::Quaterniond array_to_eigen_quat(const std::array<float, 4> &q)
{
	return Eigen::Quaterniond(q[0], q[1], q[2], q[3]);
}

double quaternion_get_yaw(const Eigen::Quaterniond &q)
{
	const double &q0 = q.w();
	const double &q1 = q.x();
	const double &q2 = q.y();
	const double &q3 = q.z();

	return std::atan2(2. * (q0 * q3 + q1 * q2), 1. - 2. * (q2 * q2 + q3 * q3));
}

} // namespace quaternion

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

#endif // !UTILS_HPP
