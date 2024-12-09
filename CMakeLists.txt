cmake_minimum_required(VERSION 3.20)

option(KN_BUILD_TESTS "Build test programs" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
        spark
        LANGUAGES CXX)

include(cmake/external.cmake)
include(cmake/random.cmake)

add_library(spark STATIC
        "src/random/seed.cpp"
        "src/particle/species.cpp"
        "src/particle/pusher.cpp"
        "src/particle/boundary.cpp"
        "src/spatial/grid.cpp"
        "src/interpolate/weight.cpp"
        "src/interpolate/field.cpp"
        "src/electromagnetics/struct_poisson.cpp"
        "src/collisions/mcc.cpp"
        "src/collisions/scattering.cpp"
        src/electromagnetics/thomas_poisson.cpp
        src/electromagnetics/util.cpp
        include/spark/core/matrix.h
)

random_src(spark)

add_library(spark::spark ALIAS spark)

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

if (KN_BUILD_TESTS)

    # CPMAddPackage("gh:catchorg/Catch2@3.5.0")
    # list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

    # enable_testing()
    # include(CTest)
    # include(Catch)

    # add_subdirectory(tests)

    # target_link_libraries(kn_tests PRIVATE
    #   spark
    #   Catch2::Catch2WithMain
    # )

endif ()
