#pragma once

#include <stdint.h>

namespace kn::random {

// Generates a pseudo-random seed based on std::time
uint64_t gen_seed();

// Initialize pseudo-random generator
void initialize(uint64_t seed);

// Uniform random uint64
uint64_t uniform_u64();

// Uniform random double in the range [0,1]
double uniform();

// Normally distributed random double with mean of 0 and standard deviation of 1
double normal();

// Normally distributed random double with specified mean and standard deviation
inline double normal(double mean, double std) {
    return std * normal() + mean;
}

// Uniform random number of type T in the range [0, vmax]
template <typename T>
inline T uniform(T vmax) {
    return static_cast<T>(kn::random::uniform() * (double)vmax);
}

// Uniform random number of type T in the range [vmin, vmax]
template <typename T>
inline T uniform(T vmin, T vmax) {
    return static_cast<T>(kn::random::uniform() * (double)(vmax - vmin)) + vmin;
}

}  // namespace kn::random
