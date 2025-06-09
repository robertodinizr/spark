#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "spark/collisions/target.h"
#include "spark/spatial/grid.h"

using namespace spark::collisions;
using namespace spark::core;
using namespace spark::spatial;


TEST_CASE("StaticUniformTarget Tests") {
    const double expected_density = 1e18;
    const double expected_temp = 300.0;

    StaticUniformTarget<2, 3> target(expected_density, expected_temp);

    SECTION("Properties are correctly returned") {
        REQUIRE(target.dens_max() == Catch::Approx(expected_density));
        REQUIRE(target.temperature() == Catch::Approx(expected_temp));
    }

    SECTION("Density is uniform everywhere") {
        REQUIRE(target.dens_at({0.0, 0.0}) == Catch::Approx(expected_density));
        REQUIRE(target.dens_at({-10.0, 50.0}) == Catch::Approx(expected_density));
        REQUIRE(target.dens_at({1e5, -1e5}) == Catch::Approx(expected_density));
    }
}


TEST_CASE("StaticFieldTarget Tests") {
    const ULongVec<2> n = {11, 11};
    const Vec<2> l = {10.0, 10.0};
    UniformGrid<2> density_grid(l, n);
    
    double min_dens = 1.0e16;
    double max_dens = 2.0e16;

    for (size_t i = 0; i < n.x; ++i) {
        for (size_t j = 0; j < n.y; ++j) {
            double fraction = static_cast<double>(i) / (n.x - 1);
            density_grid.data()(i, j) = min_dens + (max_dens - min_dens) * fraction;
        }
    }
    
    const double expected_temp = 500.0;
    StaticFieldTarget<2, 3> target(density_grid, expected_temp);

    SECTION("Properties are correctly returned") {
        REQUIRE(target.dens_max() == Catch::Approx(max_dens));
        REQUIRE(target.temperature() == Catch::Approx(expected_temp));
    }

    SECTION("Density is correctly interpolated at various positions inside grid") {
        // (i=5)
        // dens = 1e16 + (1e16) * (5 / 10) = 1.5e16
        REQUIRE(target.dens_at({5.0, 5.0}) == Catch::Approx(1.5e16));

        // (i=2.5)
        // dens = 1e16 + (1e16) * (2.5 / 10) = 1.25e16
        REQUIRE(target.dens_at({2.5, 5.0}) == Catch::Approx(1.25e16));
    }
}