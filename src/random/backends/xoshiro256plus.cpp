#ifdef KN_RANDOM_USE_XOSHIRO256PLUS

#include "random/random.h"

namespace _xoshiro256plus
{
    #include "external/xoshiro256plus.c"   
}

namespace _splitmix64
{
    #include "external/splitmix64.c"
}

namespace kn::random
{
    void initialize(uint64_t seed)
    {
        // splitmix64 is used to seed xoshiro256+ as was advised by its authors
        _splitmix64::x = seed;

        _xoshiro256plus::s[0] = _splitmix64::next();
        _xoshiro256plus::s[1] = _splitmix64::next();
        _xoshiro256plus::s[2] = _splitmix64::next();
        _xoshiro256plus::s[3] = _splitmix64::next();
    }

    uint64_t uniform_u64()
    {
        return _xoshiro256plus::next();
    }

    #define RANDOM_USE_SIMPLE_DOUBLE_CONVERSION
    double uniform()
    {   
        #ifdef RANDOM_USE_SIMPLE_DOUBLE_CONVERSION
        return (_xoshiro256plus::next() >> 11) * 0x1.0p-53;
        #else
        const union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | _xoshiro256plus::next() >> 12 };
        return u.d - 1.0;
        #endif 
    }
}

#endif