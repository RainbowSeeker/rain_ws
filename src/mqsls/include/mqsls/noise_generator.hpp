#ifndef NOISE_GENERATOR_HPP
#define NOISE_GENERATOR_HPP
#include <random>

class NoiseGenerator {
public:
    NoiseGenerator(double mean = 0.0, double stddev = 1.0)
        : gen(rd()), dist(mean, stddev) {}

    double generate() {
        return dist(gen);
    }

private:
    std::random_device rd;
    std::mt19937 gen;
    std::normal_distribution<> dist;
};

#endif // !NOISE_GENERATOR_HPP

