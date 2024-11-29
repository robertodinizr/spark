# spark ⚡

The spark (**s**tatistical **par**ticle **k**inetics) toolbox is a C++ library used for implementation of kinetic simulations with Monte Carlo-based algorithms such as the Particle-In-Cell (PIC) or Direct Simulation Monte Carlo (DSMC) methods. It aims to provide the basic functionalities for essential parts of the simulation such as integrating particle movement, solving electromagnetic field equations, and performing collisions.

> [!WARNING]
> This library is still in a very early stage of development and will suffer breaking changes.

## Building

For building, you will need CMake and a compiler that supports at least C++20. All dependencies of the toolbox are fetched automatically using [CPM](https://github.com/cpm-cmake/CPM.cmake).

To integrate the library with your project you have a couple of options. The first is to use, CMake's FetchContent as

```cmake
# In your CMakeLists.txt file

FetchContent_Declare(
  spark
  GIT_REPOSITORY https://github.com/lase-unb/spark.git
  GIT_TAG        main #or a specific commit/branch.
)

FetchContent_MakeAvailable(spark)

target_link_libraries(your_target PRIVATE spark::spark)
```

Or if you add it as a subdirectory (as in the case of a git submodule) you can simply do 

```cmake
# In your CMakeLists.txt file

add_subdirectory(spark)

target_link_libraries(your_target PRIVATE spark::spark)
```

## Dependencies 

The dependencies used by spark which are automatically downloaded and linked by CMake are 

- [HYPRE](https://github.com/hypre-space/hypre)
- [parallel-hashmap](https://github.com/greg7mdp/parallel-hashmap)
- [Eigen](https://gitlab.com/libeigen/eigen)