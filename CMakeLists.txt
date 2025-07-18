cmake_minimum_required(VERSION 3.20)

option(SPARK_BUILD_TESTS "Build test programs" ON)

option(SPARK_ENABLE_LOG_DEBUG "Log debug messages" OFF)
option(SPARK_ENABLE_LOG_INFO "Log info messages" OFF)
option(SPARK_ENABLE_LOG_WARN "Log warning messages" ON)
option(SPARK_ENABLE_LOG_ERROR "Log error messages" ON)
option(SPARK_ENABLE_LOG_ALL "Log all messages" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
        spark
        LANGUAGES CXX)

include(cmake/external.cmake)
include(cmake/random.cmake)

add_library(spark STATIC
        src/random/seed.cpp
        src/particle/pusher.cpp
        src/particle/boundary.cpp
        src/spatial/grid.cpp
        src/interpolate/weight.cpp
        src/interpolate/field.cpp
        src/em/struct_poisson.cpp
        src/collisions/mcc.cpp
        src/collisions/scattering.cpp
        src/em/thomas_poisson.cpp
        src/em/util.cpp
        src/em/electric_field.cpp
        src/particle/tiled_boundary.cpp
        src/em/cylindrical_poisson.cpp
)

if (SPARK_ENABLE_LOG_DEBUG OR SPARK_LOG_ALL)
    list(APPEND SPARK_LOG_DEFS "SPARK_ENABLE_LOG_DEBUG")
endif ()

if (SPARK_ENABLE_LOG_INFO OR SPARK_LOG_ALL)
    list(APPEND SPARK_LOG_DEFS "SPARK_ENABLE_LOG_INFO")
endif ()

if (SPARK_ENABLE_LOG_WARN OR SPARK_LOG_ALL)
    list(APPEND SPARK_LOG_DEFS "SPARK_ENABLE_LOG_WARN")
endif ()

if (SPARK_ENABLE_LOG_ERROR OR SPARK_LOG_ALL)
    list(APPEND SPARK_LOG_DEFS "SPARK_ENABLE_LOG_ERROR")
endif ()

random_src(spark)

add_library(spark::spark ALIAS spark)

target_compile_definitions(spark PRIVATE
        ${SPARK_LOG_DEFS}
)

target_include_directories(spark
        PUBLIC
        "include"
        PRIVATE
        "src"
)

target_link_libraries(spark
        PRIVATE
        phmap
        HYPRE
)

enable_testing()

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0
)

FetchContent_MakeAvailable(Catch2)

include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)

add_executable(tests 
               tests/spatial/test_grid.cpp
               tests/core/test_vec.cpp
               tests/core/test_matrix.cpp
               tests/random/test_random.cpp
               tests/particle/test_species.cpp
               
               tests/particle/test_boundary.cpp
               tests/particle/test_emitter.cpp
               tests/em/test_poisson.cpp
               tests/em/test_electric_field.cpp
               tests/interpolate/test_weight.cpp
               tests/interpolate/test_field.cpp
               tests/collisions/test_target.cpp
               tests/collisions/test_scattering.cpp
               tests/collisions/test_reactions.cpp)
target_link_libraries(tests PRIVATE spark::spark Catch2::Catch2WithMain)

catch_discover_tests(tests)
