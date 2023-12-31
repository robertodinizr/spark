cmake_minimum_required(VERSION 3.20)

option(SCI_BUILD_TESTS "Build test programs" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
  sci
  LANGUAGES CXX)

include(cmake/external_libs.cmake)

add_library(sci INTERFACE
"src/sci.h"
"src/interpolate/interpolate.h"
)

target_include_directories(sci INTERFACE
"src"
)

if(SCI_BUILD_TESTS)

enable_testing()
include(CTest)
include(Catch)

add_subdirectory(tests)

target_link_libraries(sci_tests PRIVATE 
  sci 
  Catch2::Catch2WithMain
)

endif()
