#ifdef KN_RANDOM_USE_STD_MT19937

#include "spark/random/random.h"

#include <random>

namespace _std_mt19937_64
{
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::uniform_int_distribution<uint64_t> uniform_u64;
    std::normal_distribution<double> normal;
}

namespace spark::random
{
    void initialize(uint64_t seed)
    {
        _std_mt19937_64::gen.seed(seed);
    }

    uint64_t uniform_u64()
    {
        return _std_mt19937_64::uniform_u64(_std_mt19937_64::gen);
    }

    double uniform()
    {   
        return _std_mt19937_64::uniform(_std_mt19937_64::gen);
    }

    double normal()
    {
        return _std_mt19937_64::normal(_std_mt19937_64::gen);
    }
}

#endif