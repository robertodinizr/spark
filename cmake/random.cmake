option(KN_RANDOM_USE_XOSHIRO256PLUS "Use xoroshiro256+ for pseudo-random generation" ON)
option(KN_RANDOM_USE_SPLITMIX64 "Use splitmix64 for pseudo-random generation" OFF)
option(KN_RANDOM_USE_STD_MT19937 "Use std marsenne twister (MT19937) for pseudo-random generation" OFF)

option(KN_RANDOM_SEED_USE_TIME_SPLITMIX64 "Use current time and splitmix64 to generate random seed" ON)
option(KN_RANDOM_SEED_USE_STD "Use std::random_device to generate random seed" OFF)

function(random_src tgt)
    
    # PRNG source code 
    if(KN_RANDOM_USE_XOSHIRO256PLUS)

        target_sources(${tgt} PRIVATE 
        "src/random/backends/xoshiro256plus.cpp"
        "src/random/backends/ziggurat.cpp"
        )
    
        target_compile_definitions(${tgt} PRIVATE KN_RANDOM_USE_XOSHIRO256PLUS)

    elseif(KN_RANDOM_USE_SPLITMIX64)
        
        target_sources(${tgt} PRIVATE 
        "src/random/backends/splitmix64.cpp"
        "src/random/backends/ziggurat.cpp"
        )

        target_compile_definitions(${tgt} PRIVATE KN_RANDOM_USE_SPLITMIX64)

    elseif(KN_RANDOM_USE_STD_MT19937)

        target_sources(${tgt} PRIVATE 
        "src/random/backends/std_mt19937.cpp"
        )

        target_compile_definitions(${tgt} PRIVATE KN_RANDOM_USE_STD_MT19937)    

    endif()

    if(KN_RANDOM_SEED_USE_TIME_SPLITMIX64)
        target_compile_definitions(${tgt} PRIVATE KN_RANDOM_SEED_USE_TIME_SPLITMIX64)   
    elseif(KN_RANDOM_SEED_USE_STD)
        target_compile_definitions(${tgt} PRIVATE KN_RANDOM_SEED_USE_STD)   
    endif()

endfunction(random_src)


