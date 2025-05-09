#pragma once

#include <Eigen/Eigen>
#include <Eigen/Geometry>

class TrajectoryGenerator
{
public:
    struct traj_out
    {
        Eigen::Vector3d position;
        Eigen::Vector3d velocity;
        Eigen::Vector3d acceleration;
    };
    
    virtual void update(const uint64_t dt, traj_out &out) = 0;
};

class PointTrajectoryGenerator : public TrajectoryGenerator
{
public:
    PointTrajectoryGenerator(const Eigen::Vector3d &position)
        : _position(position)
    {
    }

    void update(const uint64_t dt, traj_out &out) override
    {
        out.position = _position;
        out.velocity = Eigen::Vector3d::Zero();
        out.acceleration = Eigen::Vector3d::Zero();
    }

private:
    const Eigen::Vector3d _position;
};


#define usec2sec(x) ((double)(x) / 1e6)

class LineTrajectoryGenerator : public TrajectoryGenerator
{
public:
    LineTrajectoryGenerator(const Eigen::Vector3d &start, const Eigen::Vector3d &end, const double average_speed)
        : _duration(Eigen::Vector3d(end - start).norm() / average_speed), _avg_speed((end - start) / _duration), _current_position(start)
    {
    }

    void update(const uint64_t dt, traj_out &out) override
    {
        if (_passed_time <= _duration)
        {
            _passed_time += usec2sec(dt);
            _current_position += _avg_speed * usec2sec(dt);
            _current_velocity = _avg_speed;
        }
        else
        {
            _current_velocity = Eigen::Vector3d::Zero();
        }
        out.position = _current_position;
        out.velocity = _current_velocity;
        out.acceleration = Eigen::Vector3d::Zero();
    }

private:
    const double _duration;
    const Eigen::Vector3d _avg_speed;

    Eigen::Vector3d _current_position;
    Eigen::Vector3d _current_velocity;

    double _passed_time = 0.0;
};

class CircleTrajectoryGenerator : public TrajectoryGenerator
{
public:
    CircleTrajectoryGenerator(const Eigen::Vector3d &center, const double radius, const double omega)
        : _center(center), _radius(radius), _omega(omega)
    {
    }

    void update(const uint64_t dt, traj_out &out) override
    {
        _passed_time += usec2sec(dt);
        const double theta = _omega * _passed_time;
        out.position = _center + Eigen::Vector3d(_radius * cos(theta), _radius * sin(theta), 0);
        out.velocity = Eigen::Vector3d(-_radius * _omega * sin(theta), _radius * _omega * cos(theta), 0);
        // out.acceleration = Eigen::Vector3d::Zero();
        out.acceleration = Eigen::Vector3d(-_radius * _omega * _omega * cos(theta), -_radius * _omega * _omega * sin(theta), 0);
    }

private:
    const Eigen::Vector3d _center;
    const double _radius;
    const double _omega;

    double _passed_time = 0.0;
};

class RectangleTrajectoryGenerator : public TrajectoryGenerator
{
    enum class Phase
    {
        L_TO_R,
        U_TO_D,
        R_TO_L,
        D_TO_U
    };
public:
    RectangleTrajectoryGenerator(const Eigen::Vector3d &start, const Eigen::Vector3d &end, const double average_speed)
        : _start(start), _end(end), _average_speed(average_speed), _current_position(start)
    {
    }

    void update(const uint64_t dt, traj_out &out) override
    {
        switch (_phase)
        {
        case Phase::L_TO_R:
            _current_position[0] += _average_speed * usec2sec(dt);
            if (_current_position[0] >= _end[0])
            {
                _current_position[0] = _end[0];
                _phase = Phase::U_TO_D;
            }
            _current_velocity = {_average_speed, 0, 0};
            break;
        case Phase::U_TO_D:
            _current_position[1] += _average_speed * usec2sec(dt);
            if (_current_position[1] >= _end[1])
            {
                _current_position[1] = _end[1];
                _phase = Phase::R_TO_L;
            }
            _current_velocity = {0, _average_speed, 0};
            break;
        case Phase::R_TO_L:
            _current_position[0] -= _average_speed * usec2sec(dt);
            if (_current_position[0] <= _start[0])
            {
                _current_position[0] = _start[0];
                _phase = Phase::D_TO_U;
            }
            _current_velocity = {-_average_speed, 0, 0};
            break;
        case Phase::D_TO_U:
            _current_position[1] -= _average_speed * usec2sec(dt);
            if (_current_position[1] <= _start[1])
            {
                _current_position[1] = _start[1];
                _phase = Phase::L_TO_R;
            }
            _current_velocity = {0, -_average_speed, 0};
            break;
        }
        out.position = _current_position;
        out.velocity = _current_velocity;
        out.acceleration = Eigen::Vector3d::Zero();
    }

private:
    const Eigen::Vector3d _start;
    const Eigen::Vector3d _end;
    const double _average_speed;

    Phase _phase = Phase::L_TO_R;
    Eigen::Vector3d _current_position;
    Eigen::Vector3d _current_velocity;
};

class LissajousTrajectoryGenerator : public TrajectoryGenerator
{
public:
    /**
     * @brief Construct a new Lissajous Trajectory Generator object
     * @param center Center of the Lissajous curve
     * @param a Amplitude of x-axis
     * @param b Amplitude of y-axis
     * @param omega_x Angular frequency of x-axis
     * @param omega_y Angular frequency of y-axis
     * @param phi Phase difference
     */
    LissajousTrajectoryGenerator(const Eigen::Vector3d &center, const double a, const double b, const double omega_x, const double omega_y, const double phi)
        : _center(center), _a(a), _b(b), _omega_x(omega_x), _omega_y(omega_y), _phi(phi)
    {
    }

    void update(const uint64_t dt, traj_out &out) override
    {
        _passed_time += usec2sec(dt);
        const double x = _a * sin(_omega_x * _passed_time + _phi);
        const double y = _b * sin(_omega_y * _passed_time);
        out.position = _center + Eigen::Vector3d(x, y, 0);
        out.velocity = Eigen::Vector3d(_a * _omega_x * cos(_omega_x * _passed_time + _phi), _b * _omega_y * cos(_omega_y * _passed_time), 0);
        out.acceleration = Eigen::Vector3d(-_a * _omega_x * _omega_x * sin(_omega_x * _passed_time + _phi), -_b * _omega_y * _omega_y * sin(_omega_y * _passed_time), 0);
    }

private:
    const Eigen::Vector3d _center;
    const double _a, _b;
    const double _omega_x, _omega_y;
    const double _phi;

    double _passed_time = 0.0;
};

static std::shared_ptr<TrajectoryGenerator> make_trajectory_generator(const std::string &traj_type, const std::string &test_mode)
{
    #define deg2rad(x)  ((x) * M_PI / 180)

    double HGT_DEFAULT = -0.3;
    if (test_mode == "sil" || test_mode == "hil") {
        HGT_DEFAULT = -1;
    }

    if (traj_type == "point") {
        return std::make_shared<PointTrajectoryGenerator>(Eigen::Vector3d(0, 0, HGT_DEFAULT));
    } else if (traj_type == "line") {
        return std::make_shared<LineTrajectoryGenerator>(Eigen::Vector3d(2, 0, HGT_DEFAULT), Eigen::Vector3d(-2, 0, HGT_DEFAULT), 0.25);
    } else if (traj_type == "circle") {
        return std::make_shared<CircleTrajectoryGenerator>(Eigen::Vector3d(0, 0, HGT_DEFAULT), 10, deg2rad(12));
    } else if (traj_type == "rectangle") {
        return std::make_shared<RectangleTrajectoryGenerator>(Eigen::Vector3d(0, 0, HGT_DEFAULT), Eigen::Vector3d(20, 20, HGT_DEFAULT), 2);
    } else if (traj_type == "lissajous") {
        return std::make_shared<LissajousTrajectoryGenerator>(Eigen::Vector3d(0, 0, HGT_DEFAULT), 2, 1, deg2rad(12), deg2rad(24), deg2rad(90));
    } else {
        std::cerr << "Invalid trajectory type: " << traj_type << std::endl;
        return nullptr;
    }
}