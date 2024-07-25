#include "kn/random/random.h"


#ifdef KN_RANDOM_SEED_USE_TIME_SPLITMIX64

#include <ctime>
#include <mutex>

namespace _seed
{
    #include "backends/external/splitmix64.c"    
}

namespace
{
    std::once_flag once;
}

namespace kn::random
{
    uint64_t gen_seed()
    {
        std::call_once(once, [&](){
            _seed::x = std::time(NULL);
        });

        return _seed::next();
    }
}

#endif

#ifdef KN_RANDOM_SEED_USE_STD

#include <random>

namespace _seed
{
    std::random_device device;
    std::uniform_int_distribution<uint64_t> uniform_u64;
}

namespace kn::random
{
    uint64_t gen_seed()
    {
        return _seed::uniform_u64(_seed::device);
    }
}


#endif