cmake_minimum_required(VERSION 3.20)

option(KN_BUILD_TESTS "Build test programs" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
  kn
  LANGUAGES CXX)

include(cmake/external.cmake)
include(cmake/random.cmake)

add_library(kn STATIC
  "src/random/seed.cpp"
  "src/particle/species.cpp"
  "src/particle/pusher.cpp"
  "src/particle/boundary.cpp"
  "src/spatial/grid.cpp"
  "src/interpolate/weight.cpp"
  "src/interpolate/field.cpp"
  "src/electromagnetics/poisson.cpp"
  "src/collisions/mcc.cpp"
  "src/collisions/reaction.cpp"
  "src/collisions/scattering.cpp"
)

random_src(kn)

add_library(kn::kn ALIAS kn)

target_include_directories(kn 
PUBLIC
  "include"
PRIVATE
  "src"
)

target_link_libraries(kn PRIVATE Eigen)

if(KN_BUILD_TESTS)

# CPMAddPackage("gh:catchorg/Catch2@3.5.0")
# list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

# enable_testing()
# include(CTest)
# include(Catch)

# add_subdirectory(tests)

# target_link_libraries(kn_tests PRIVATE 
#   kn 
#   Catch2::Catch2WithMain
# )

endif()
