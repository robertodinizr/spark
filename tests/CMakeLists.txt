add_executable(kn_tests test_grid.cpp)

target_include_directories(kn_tests PRIVATE
    ${Catch2_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(kn_tests PRIVATE spark_lib)

catch_discover_tests(kn_tests)