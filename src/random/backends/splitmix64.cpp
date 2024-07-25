#ifdef KN_RANDOM_USE_SPLITMIX64

#include "random/random.h"

namespace _splitmix64
{
    #include "external/splitmix64.c"
}

namespace kn::random
{
    void initialize(uint64_t seed)
    {
        _splitmix64::x = seed;
    }

    uint64_t uniform_u64()
    {
        return _splitmix64::next();
    }

    #define RANDOM_USE_SIMPLE_DOUBLE_CONVERSION
    double uniform()
    {   
        #ifdef RANDOM_USE_SIMPLE_DOUBLE_CONVERSION
        return (_splitmix64::next() >> 11) * 0x1.0p-53;
        #else
        const union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | _splitmix64::next() >> 12 };
        return u.d - 1.0;
        #endif 
    }

}

#endif