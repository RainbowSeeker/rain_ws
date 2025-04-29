#pragma once

struct RefTranslation
{
    double delta_xyz[3];

    double delta_x() const
    {
        return delta_xyz[0];
    }

    double delta_y() const
    {
        return delta_xyz[1];
    }

    double delta_z() const
    {
        return delta_xyz[2];
    }

    const RefTranslation operator+(const RefTranslation &other) const
    {
        return {delta_xyz[0] + other.delta_xyz[0],
                delta_xyz[1] + other.delta_xyz[1],
                delta_xyz[2] + other.delta_xyz[2]};
    }
    const RefTranslation operator-(const RefTranslation &other) const
    {
        return {delta_xyz[0] - other.delta_xyz[0],
                delta_xyz[1] - other.delta_xyz[1],
                delta_xyz[2] - other.delta_xyz[2]};
    }
    const RefTranslation operator/(const double &other) const
    {
        return {delta_xyz[0] / other,
                delta_xyz[1] / other,
                delta_xyz[2] / other};
    }
};

class RTKComponent
{
public:
    RCLCPP_SMART_PTR_DEFINITIONS(RTKComponent)

    /*
    * @brief Handle the received message
    * @param bytes The received bytes
    * @param result The parsed result
    * @return 0 if continue, 1 if finished, -1 if error
    */
    virtual int handle_recv_message(const std::string &bytes, RefTranslation &result) = 0;

    virtual ~RTKComponent() = default;
};

