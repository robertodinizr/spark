add_executable(spark_tests STATIC
               test_grid.cpp
               test_vec.cpp
               test_matrix.cpp
               test_random.cpp)

target_include_directories(spark_tests PRIVATE
    ${Catch2_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(spark_tests PRIVATE spark_lib)

catch_discover_tests(spark_tests)