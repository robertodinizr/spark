#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "spark/interpolate/field.h"
#include "spark/particle/species.h"
#include "spark/spatial/grid.h"

using namespace spark;
using namespace spark::core;
using namespace spark::spatial;
using namespace spark::particle;
using namespace spark::interpolate;

using TestSpecies1D3V = ChargedSpecies<1, 3>;
using TestSpecies2D3V = ChargedSpecies<2, 3>;

TEST_CASE("Field to particles interpolation (1D)") {
    UniformGrid<1> grid({4.0}, {5});
    for(size_t i = 0; i < grid.n().x; ++i) {
        grid.data()[i] = 10.0 * i;
    }
    
    TestSpecies1D3V species(1.0, 1.0);
    TMatrix<double, 1> particle_field;

    SECTION("Single particle") {
        species.add(1, [&](auto& v, auto& x){ x = {1.5}; });
        field_at_particles(grid, species, particle_field);

        // Normalized position: xp_norm = 1.5
        // il = 1
        // phi_p = phi[1]*(1 - 0.5) + phi[2]*(0.5) = 10 * 0.5 + 20 * 0.5 = 5 + 10 = 15.0
        REQUIRE(particle_field[0] == Catch::Approx(15.0));
    }

    SECTION("Particle over node") {
        species.add(1, [&](auto& v, auto& x){ x = {2.0}; });
        field_at_particles(grid, species, particle_field);

        // phi_p = phi[2] = 20.0.
        REQUIRE(particle_field[0] == Catch::Approx(20.0));
    }
}

TEST_CASE("Field to particles interpolation (2D Cartesian)") {
    UniformGrid<2> grid({4.0, 4.0}, {5, 5});
    for(size_t i = 0; i < grid.n().x; ++i) {
        for(size_t j = 0; j < grid.n().y; ++j) {
            grid.data()(i,j) = static_cast<double>(i + 10*j);
        }
    }
    
    TestSpecies2D3V species(1.0, 1.0);
    TMatrix<double, 1> particle_field;

    SECTION("Single particle") {
        species.add(1, [&](auto& v, auto& x){ x = {1.2, 2.7}; });
        field_at_particles(grid, species, particle_field);
        
        // x_local = 0.2, y_local = 0.7
        // phi(1,2)=21, phi(2,2)=22, phi(1,3)=31, phi(2,3)=32
        // w(1,2) = (1-0.2)*(1-0.7) = 0.24
        // w(2,2) = 0.2*(1-0.7) = 0.06
        // w(1,3) = (1-0.2)*0.7 = 0.56
        // w(2,3) = 0.2*0.7 = 0.14
        // phi_p = 21*0.24 + 22*0.06 + 31*0.56 + 32*0.14 = 5.04 + 1.32 + 17.36 + 4.48 = 28.2
        REQUIRE(particle_field[0] == Catch::Approx(28.2));
    }
}

TEST_CASE("Field to particles interpolation (2D Cylindrical)") {
    UniformGrid<2> grid({4.0, 4.0}, {5, 5});
    for(size_t i = 0; i < grid.n().x; ++i) {
        for(size_t j = 0; j < grid.n().y; ++j) {
            grid.data()(i,j) = static_cast<double>(i + 10*j);
        }
    }
    
    TestSpecies2D3V species(1.0, 1.0);
    TMatrix<double, 1> particle_field;

    SECTION("Single particle") {
        species.add(1, [&](auto& v, auto& x){ x = {1.2, 2.7}; });
        field_at_particles_cylindrical(grid, species, particle_field);
        
        // a1=0.2736, a2=0.0684, a3=0.1316, a4=0.5264
        // phi(1,2)=21, phi(2,2)=22, phi(1,3)=31, phi(2,3)=32
        // phi_p = phi(1,2)*a1 + phi(2,2)*a2 + phi(2,3)*a3 + phi(1,3)*a4 = 21*0.2736 + 22*0.0684 + 32*0.1316 + 31*0.5264 = 5.7456 + 1.5048 + 4.2112 + 16.3184 = 27.78
        REQUIRE(particle_field[0] == Catch::Approx(27.78));
    }
    
    SECTION("Particle at r = 0") {
        species.add(1, [&](auto& v, auto& x){ x = {2.5, 0.0}; });
        field_at_particles_cylindrical(grid, species, particle_field);
        
        // z_local = 0.5, r_local = 0.0
        // a1=0.5, a2=0.5
        // phi_p = phi(2,0)*a1 + phi(3,0)*a2 = 2*0.5 + 3*0.5 = 2.5
        REQUIRE(particle_field[0] == Catch::Approx(2.5));
    }

    SECTION("Particle with r < 0") {
        species.add(1, [&](auto& v, auto& x){ x = {2.5, -0.1}; });
        field_at_particles_cylindrical(grid, species, particle_field);
        REQUIRE(particle_field[0] == Catch::Approx(0.0));
    }
}