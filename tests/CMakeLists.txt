add_executable(spark_tests STATIC
               spatial/test_grid.cpp
               core/test_vec.cpp
               core/test_matrix.cpp
               random/test_random.cpp
               particle/test_species.cpp
               particle/test_pusher.cpp
               particle/test_boundary.cpp
               em/test_poisson.cpp
               em/test_electric_field.cpp
               interpolate/test_weight.cpp
               interpolate/test_field.cpp
               collisions/test_target.cpp
               collisions/test_scattering.cpp
               collisions/test_reactions.cpp
               collisions/test_mcc.cpp)

target_include_directories(spark_tests PRIVATE
    ${Catch2_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(spark_tests PRIVATE spark_lib)

catch_discover_tests(spark_tests)