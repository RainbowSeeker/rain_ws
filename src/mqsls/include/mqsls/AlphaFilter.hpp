
/**
 * @file AlphaFilter.hpp
 *
 * @brief First order "alpha" IIR digital filter also known as leaky integrator or forgetting average.
 *
 * @author Mathieu Bresciani <brescianimathieu@gmail.com>
 * @author Matthias Grob <maetugr@gmail.com>
 */

#pragma once

#include <cfloat>
#include <cmath>

template <typename T>
class AlphaFilter
{
public:
	AlphaFilter() = default;
	explicit AlphaFilter(double alpha) : _alpha(alpha) {}
	explicit AlphaFilter(double alpha, const T &initial_state) : _alpha(alpha), _filter_state(initial_state) {}

	~AlphaFilter() = default;

	/**
	 * Set filter parameters for time abstraction
	 *
	 * Both parameters have to be provided in the same units.
	 *
	 * @param sample_interval interval between two samples
	 * @param time_constant filter time constant determining convergence
	 */
	void setParameters(double sample_interval, double time_constant)
	{
		const double denominator = time_constant + sample_interval;

		if (denominator > DBL_EPSILON) {
			setAlpha(sample_interval / denominator);
		}
	}

	bool setCutoffFreq(double sample_freq, double cutoff_freq)
	{
		if ((sample_freq <= 0.0) || (cutoff_freq <= 0.0) || (cutoff_freq >= sample_freq / 2.0)
		    || !std::isfinite(sample_freq) || !std::isfinite(cutoff_freq)) {

			// Invalid parameters
			return false;
		}

		setParameters(1.0 / sample_freq, 1.0 / (2.0 * M_PI * cutoff_freq));
		_cutoff_freq = cutoff_freq;
		return true;
	}

	/**
	 * Set filter parameter alpha directly without time abstraction
	 *
	 * @param alpha [0,1] filter weight for the previous state. High value - long time constant.
	 */
	void setAlpha(double alpha) { _alpha = alpha; }

	/**
	 * Set filter state to an initial value
	 *
	 * @param sample new initial value
	 */
	void reset(const T &sample) { _filter_state = sample; }

	/**
	 * Add a new raw value to the filter
	 *
	 * @return retrieve the filtered result
	 */
	const T &update(const T &sample)
	{
		_filter_state = updateCalculation(sample);
		return _filter_state;
	}

	const T &getState() const { return _filter_state; }
	double getCutoffFreq() const { return _cutoff_freq; }

protected:
	T updateCalculation(const T &sample) { return (1.0 - _alpha) * _filter_state + _alpha * sample; }

	double _cutoff_freq{0.0};
	double _alpha{0.0};
	T _filter_state{};
};
