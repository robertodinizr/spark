 add_executable(kn_tests test_grid.cpp)
target_include_directories(kn_tests PRIVATE ${Catch2_SOURCE_DIR}/src)
 catch_discover_tests(kn_tests)