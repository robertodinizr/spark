#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "spark/interpolate/weight.h"
#include "spark/particle/species.h"
#include "spark/spatial/grid.h"

using namespace spark;
using namespace spark::core;
using namespace spark::spatial;
using namespace spark::particle;
using namespace spark::interpolate;

using TestSpecies1D3V = ChargedSpecies<1, 3>;
using TestSpecies2D3V = ChargedSpecies<2, 3>;

TEST_CASE("Weighting particles to grid (1D)") {
    UniformGrid<1> grid({4.0}, {5});
    TestSpecies1D3V species(1.0, 1.0);

    SECTION("Single particle") {
        species.add(1, [&](auto& v, auto& x){ x = {1.5}; });
        weight_to_grid(species, grid);

        // Normalized position: xp_dx = x_p / dx = 1.5 / 1.0 = 1.5
        // Left node: il = floor(1.5) = 1
        // Weighting left node: w_left = il + 1 - xp_dx = 1.0 + 1.0 - 1.5 = 0.5
        // Weighting right node: w_right = xp_dx - il = 1.5 - 1.0 = 0.5
        auto& data = grid.data();
        REQUIRE(data[1] == Catch::Approx(0.5));
        REQUIRE(data[2] == Catch::Approx(0.5));
        REQUIRE(data[0] + data[3] + data[4] == Catch::Approx(0.0));
    }

    SECTION("Multiple particles") {
        species.add(1, [&](auto& v, auto& x){ x = {1.5}; });
        species.add(1, [&](auto& v, auto& x){ x = {1.8}; });
        weight_to_grid(species, grid);

        // Node 1: 0.5 particle_1 + 0.2 particle_2 = 0.7
        // Node 2: 0.5 particle_1 + 0.8 particle_2 = 1.3
        auto& data = grid.data();
        REQUIRE(data[1] == Catch::Approx(0.7));
        REQUIRE(data[2] == Catch::Approx(1.3));
    }
}

TEST_CASE("Weighting particles to grid (2D Cartesian)") {
    UniformGrid<2> grid({4.0, 4.0}, {5, 5});
    TestSpecies2D3V species(1.0, 1.0);

    SECTION("Single particle") {
        species.add(1, [&](auto& v, auto& x){ x = {1.2, 2.7}; });
        weight_to_grid(species, grid);

        // Normalized position: xp = 1.2, yp = 2.7
        // Nodes: (1,2), (2,2), (1,3), (2,3)
        // Local distance: x_local = 0.2, y_local = 0.7
        // Weighting:
        // Node (1,2) [inferior-left]: (1-x_local)*(1-y_local) = (1-0.2)*(1-0.7) = 0.8 * 0.3 = 0.24
        // Node (2,2) [inferior-right]: x_local *(1-y_local) = 0.2 *(1-0.7) = 0.2 * 0.3 = 0.06
        // Node (1,3) [superior-left]: (1-x_local)* y_local = (1-0.2)* 0.7 = 0.8 * 0.7 = 0.56
        // Node (2,3) [superior-right]: x_local * y_local = 0.2 * 0.7 = 0.14
        // Total: 0.24 + 0.06 + 0.56 + 0.14 = 1.0
        auto& data = grid.data();
        REQUIRE(data(1,2) == Catch::Approx(0.24));
        REQUIRE(data(2,2) == Catch::Approx(0.06));
        REQUIRE(data(1,3) == Catch::Approx(0.56));
        REQUIRE(data(2,3) == Catch::Approx(0.14));
    }
}

TEST_CASE("Weighting particles to grid (2D Cylindrical)") {
    UniformGrid<2> grid({4.0, 4.0}, {5, 5});
    TestSpecies2D3V species(1.0, 1.0);
    
    SECTION("Single particle") {
        species.add(1, [&](auto& v, auto& x){ x = {1.2, 2.7}; });
        weight_to_grid_cylindrical(species, grid);
        
        // Coords (z,r): (1.2, 2.7) 
        // Normalized position: zp = 1.2, rp = 2.7
        // Nodes: (1,2), (2,2), (1,3), (2,3)
        // Local distance: z_local = 0.2, r_local = 0.7
        // rj = 2 * dr = 2.0.
        // denominator = rj + 0.5*dr = 2.5
        // f1 = (rj + 0.5*r_local*dr) / 2.5 = (2.0 + 0.5*0.7*1.0) / 2.5 = 0.94
        // f2 = (rj + 0.5*(r_local+1)*dr) / 2.5 = (2.0 + 0.5*1.7*1.0) / 2.5 = 1.14
        // Node (1,2): (1-z_local)*(1-r_local)*f2 = 0.8 * 0.3 * 1.14 = 0.2736
        // Node (2,2): z_local*(1-r_local)*f2 = 0.2 * 0.3 * 1.14 = 0.0684
        // Node (1,3): (1-z_local)*r_local*f1 = 0.8 * 0.7 * 0.94 = 0.5264
        // Node (2,3): z_local*r_local*f1 = 0.2 * 0.7 * 0.94 = 0.1316
        auto& data = grid.data();
        REQUIRE(data(1,2) == Catch::Approx(0.2736));
        REQUIRE(data(2,2) == Catch::Approx(0.0684));
        REQUIRE(data(1,3) == Catch::Approx(0.5264));
        REQUIRE(data(2,3) == Catch::Approx(0.1316));
    }

    SECTION("Particle at r = 0") {
        species.add(1, [&](auto& v, auto& x){ x = {2.5, 0.0}; });
        weight_to_grid_cylindrical(species, grid);

        // Position (z,r) = (2.5, 0.0) 
        // Normalized position: zp = 2.5, rp = 0.0
        // z_local = 0.5, r_local = 0.0
        // rj = 0.0 
        // Denominator = 0.5*dr = 0.5
        // f2 = (0 + 0.5*1.0*1.0) / 0.5 = 1.0
        // Node (2,0): (1-z_local)*(1-r_local)*f2 = 0.5 * 1.0 * 1.0 = 0.5
        // Node (3,0): z_local*(1-r_local)*f2 = 0.5 * 1.0 * 1.0 = 0.5
        auto& data = grid.data();
        REQUIRE(data(2,0) == Catch::Approx(0.5));
        REQUIRE(data(3,0) == Catch::Approx(0.5));
        REQUIRE(data(2,1) == Catch::Approx(0.0));
    }

    SECTION("Particle r < 0") {
        species.add(1, [&](auto& v, auto& x){ x = {2.5, -0.1}; });
        weight_to_grid_cylindrical(species, grid);

        double total_weight = 0.0;
        auto& data = grid.data();
        for(size_t i=0; i<grid.n().x; ++i) 
            for(size_t j=0; j<grid.n().y; ++j) 
                total_weight += data(i,j);
        REQUIRE(total_weight == Catch::Approx(0.0));
    }
}