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

#define usec2sec(x) ((double)(x) / 1e6)

class LineTrajectoryGenerator : public TrajectoryGenerator
{
public:
    LineTrajectoryGenerator(const Eigen::Vector3d &start, const Eigen::Vector3d &end, const uint64_t duration)
        : _duration(usec2sec(duration)), _avg_speed((end - start) / _duration), _current_position(start)
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
        out.acceleration = Eigen::Vector3d::Zero();//Eigen::Vector3d(-_radius * _omega * _omega * cos(theta), -_radius * _omega * _omega * sin(theta), 0);
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
    RectangleTrajectoryGenerator(const Eigen::Vector3d &start, const Eigen::Vector3d &end, const int average_speed)
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
    const int _average_speed;

    Phase _phase = Phase::L_TO_R;
    Eigen::Vector3d _current_position;
    Eigen::Vector3d _current_velocity;
};

