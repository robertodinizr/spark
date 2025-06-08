#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "spark/particle/boundary.h"
#include "spark/particle/species.h"
#include "spark/core/vec.h"
#include "spark/spatial/grid.h"

using namespace spark::particle;
using namespace spark::core;
using namespace spark::spatial;

using TestSpecies1D = ChargedSpecies<1, 3>;
using TestSpecies2D = ChargedSpecies<2, 3>;

TEST_CASE("1D Boundary Conditions") {
    const double xmin = 0.0;
    const double xmax = 10.0;

    SECTION("Absorbing Boundary") {
        TestSpecies1D species(1.0, 1.0);
        species.add(1, [&](auto& v, auto& x){ x = {5.0}; });
        species.add(1, [&](auto& v, auto& x){ x = {xmin}; });
        species.add(1, [&](auto& v, auto& x){ x = {xmax}; });
        species.add(1, [&](auto& v, auto& x){ x = {xmin - 1.0}; });
        species.add(1, [&](auto& v, auto& x){ x = {xmax + 1.0}; });

        REQUIRE(species.n() == 5);
        apply_absorbing_boundary(species, xmin, xmax);
        REQUIRE(species.n() == 3);
        REQUIRE(species.x()[0].x == 5.0);
        REQUIRE(species.x()[1].x == xmin);
        REQUIRE(species.x()[2].x == xmax);
    }

    SECTION("Symmetric Boundary") {
        TestSpecies1D species(1.0, 1.0);
        species.add(1, [&](auto& v, auto& x){ x = {xmin - 1.0}; });
        species.add(1, [&](auto& v, auto& x){ x = {xmax + 2.0}; });
        species.add(1, [&](auto& v, auto& x){ x = {5.0}; });

        apply_symmetric_boundary(species, xmin, xmax);
        REQUIRE(species.n() == 3);
        // dx = xmin - pos = 0 - (-1) = 1 
        // pos_new = xmax - dx = 10 - 1 = 9
        // Particle 2: dx = pos - xmax = 12 - 10 = 2 
        // pos_new = xmin + dx = 0 + 2 = 2
        REQUIRE(species.x()[0].x == Catch::Approx(xmax - 1.0));
        REQUIRE(species.x()[1].x == Catch::Approx(xmin + 2.0));
        REQUIRE(species.x()[2].x == Catch::Approx(5.0));
    }
}


TEST_CASE("2D Tiled Boundary") {
    const GridProp<2> grid_prop = {
        .l = {10.0, 10.0},
        .dx = {1.0, 1.0},
        .n = {11, 11}
    };
    const double dt = 1.0;

    const std::vector<TiledBoundary> boundary_defs = {
        {{1, 1}, {1, 9}, BoundaryType::Specular},
        {{9, 1}, {9, 9}, BoundaryType::Specular},
        {{1, 1}, {9, 1}, BoundaryType::Absorbing},
        {{1, 9}, {9, 9}, BoundaryType::Specular},
    };

    TiledBoundary2D boundary(grid_prop, boundary_defs, dt);
    TestSpecies2D species(1.0, 1.0);

    SECTION("Particle does not collide") {
        species.add(1, [&](auto& v, auto& x){
            v = {1.0, 1.0, 0.0};
            x = {3.0, 3.0};
        });
        
        boundary.apply(species);
        REQUIRE(species.n() == 1);
        REQUIRE(species.x()[0].x == Catch::Approx(3.0));
        REQUIRE(species.x()[0].y == Catch::Approx(3.0));
    }

    SECTION("Particle is absorbed by a wall") {
        species.add(1, [&](auto& v, auto& x){
            v = {0.0, -2.0, 0.0};
            x = {5.0, 0.5};
        });

        boundary.apply(species);
        REQUIRE(species.n() == 0);
    }
    
    SECTION("Particle reflects on a vertical wall") {
        species.add(1, [&](auto& v, auto& x){
            v = {-2.0, 1.0, 0.0};
            x = {0.5, 5.0};
        });

        boundary.apply(species);
        REQUIRE(species.n() == 1);


        // v = {-2.0, 1.0}
        // x = {0.5, 5.0}
        // dt = 1.0
        // x0 = x1 - v*dt = {0.5 - (-2.0), 5.0 - 1.0} = {2.5, 4.0}
        // x(t) = x0.x + t * v.x = 2.5 - 2.0 * t
        // y(t) = y0.y + t * v.y = 4.0 + 1.0 * t
        // t = 0.25
        // contact = {x(0.25), y(0.25)} = {2.0, 4.25}
        // x_reflected = x1 + 2*n*(contact - x1) = 0.5 + 2*1*(2.0-0.5) = 3.5
        REQUIRE(species.x()[0].x == Catch::Approx(3.5));
        REQUIRE(species.x()[0].y == Catch::Approx(5.0));
        REQUIRE(species.v()[0].x == Catch::Approx(2.0));
        REQUIRE(species.v()[0].y == Catch::Approx(1.0));
    }

    SECTION("Particle reflects on a horizontal wall") {
        species.add(1, [&](auto& v, auto& x){
            v = {1.0, 2.0, 0.0};
            x = {5.0, 9.5};
        });

        boundary.apply(species);
        REQUIRE(species.n() == 1);

        // y1 = 9.5
        // contact = {..., 9.0}
        // normal = {0.0, -1.0}
        // y_refletido = y_final + 2*n.y*(contato.y - y_final) = 9.5 + 2*(-1)*(9.0-9.5) = 8.5
        REQUIRE(species.x()[0].x == Catch::Approx(5.0));
        REQUIRE(species.x()[0].y == Catch::Approx(8.5));
        REQUIRE(species.v()[0].x == Catch::Approx(1.0));
        REQUIRE(species.v()[0].y == Catch::Approx(-2.0));
    }

    SECTION("Particle starts inside a wall and is removed") {
        species.add(1, [&](auto& v, auto& x){ v = {0.0, 0.0, 0.0}; x = {5.0, 1.0}; });
        boundary.apply(species);
        REQUIRE(species.n() == 0);
    }
}