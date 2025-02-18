#ifndef PID_HPP
#define PID_HPP

namespace mqsls
{
    
template <typename T>
class Pid
{
public:
    /*!
    * \brief Constructor, zeros out Pid values when created and
    *        initialize Pid-gains and integral term limits.
    *        Does not initialize dynamic reconfigure for PID gains
    *
    * \param p The proportional gain.
    * \param i The integral gain.
    * \param d The derivative gain.
    * \param i_max The max integral windup.
    * \param i_min The min integral windup.
    * \param antiwindup If true, antiwindup is enabled and i_max/i_min are enforced
    *
    * \throws An std::invalid_argument exception is thrown if i_min > i_max
    */
    Pid(T p = T(), T i = T(), T d = T(), T i_max = T(), T i_min = T(), bool antiwindup = false)
    {
        initPid(p, i, d, i_max, i_min, antiwindup);
    }

    /*!
    * \brief Destructor of Pid class.
    */
    ~Pid()
    {
    }

    /*!
    * \brief Zeros out Pid values and initialize Pid-gains and integral term limits
    *        Does not initialize the node's parameter interface for PID gains
    *
    * \param p The proportional gain.
    * \param i The integral gain.
    * \param d The derivative gain.
    * \param i_max The max integral windup.
    * \param i_min The min integral windup.
    * \param antiwindup If true, antiwindup is enabled and i_max/i_min are enforced
    *
    * \note New gains are not applied if i_min_ > i_max_
    */
    void initPid(T p, T i, T d, T i_max, T i_min, bool antiwindup)
    {
        p_ = p;
        i_ = i;
        d_ = d;
        i_max_ = i_max;
        i_min_ = i_min;
        antiwindup_ = antiwindup;
        reset();
    }

    /*!
    * \brief Reset the state of this PID controller
    */
    void reset()
    {
        p_error_last_ = T();
        p_error_ = T();
        i_error_ = T();
        d_error_ = T();
        cmd_ = T();
        error_dot_ = T();
    }

    /*!
    * \brief Set the PID error and compute the PID command with nonuniform time
    * step size. The derivative error is computed from the change in the error
    * and the timestep \c dt.
    *
    * \param error  Error since last call (error = target - state)
    * \param dt Change in time since last call in nanoseconds
    *
    * \returns PID command
    */
    T computeCommand(T error, uint64_t dt)
    {
        T error_dot = (error - p_error_last_) / (static_cast<double>(dt) / 1e6);
        p_error_last_ = error;
        return computeCommand(error, error_dot, dt);
    }

    /*!
    * \brief Set the PID error and compute the PID command with nonuniform
    * time step size. This also allows the user to pass in a precomputed
    * derivative error.
    *
    * \param error Error since last call (error = target - state)
    * \param error_dot d(Error)/dt since last call
    * \param dt Change in time since last call in us.
    *
    * \returns PID command
    */
    T computeCommand(T error, T error_dot, uint64_t dt)
    {
        p_error_ = error;
        i_error_ += error * (static_cast<double>(dt) / 1e6);
        d_error_ = error_dot;
        error_dot_ = error_dot;

        if constexpr (std::is_same<T, double>::value)
        {
            if (antiwindup_)
                i_error_ = std::clamp(i_error_, i_min_, i_max_);
            cmd_ = p_ * p_error_ + i_ * i_error_ + d_ * d_error_;
        }
        else if constexpr (std::is_same<T, Eigen::Vector3d>::value)
        {
            if (antiwindup_)
                i_error_ = i_error_.cwiseMax(i_min_).cwiseMin(i_max_);
            cmd_ = p_.cwiseProduct(p_error_) + i_.cwiseProduct(i_error_) + d_.cwiseProduct(d_error_);
        }
        else
        {
            static_assert(std::is_same<T, double>::value || std::is_same<T, Eigen::Vector3d>::value, "Only double and Eigen::Vector3d are supported");
        }

        return cmd_;
    }

    Pid(const Pid &source) = delete;
    Pid(Pid &&source) = delete;
    Pid &operator=(const Pid &source) = delete;
private:
    T p_, i_, d_, i_max_, i_min_;
    bool antiwindup_;

    T p_error_last_; /**< _Save position state for derivative state calculation. */
    T p_error_;      /**< Position error. */
    T i_error_;      /**< Integral of position error. */
    T d_error_;      /**< Derivative of position error. */
    T cmd_;          /**< Command to send. */
    T error_dot_;    /**< Derivative error */
};

} // namespace mqsls

#endif // !PID_HPP