#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "spark/particle/species.h"
#include "spark/core/vec.h"

using namespace spark::particle;
using namespace spark::core;

using TestSpecies1D3V = ChargedSpecies<1, 3>;
using TestSpecies2D3V = ChargedSpecies<2, 3>;
using TestSpecies3D3V = ChargedSpecies<3, 3>;


TEST_CASE("Species") {

    SECTION("Species properties") {
        const double mass = 1.0;
        const double charge = -1.6e-19;
        TestSpecies1D3V species1D3V(charge, mass);
        TestSpecies2D3V species2D3V(charge, mass);
        TestSpecies3D3V species3D3V(charge, mass);

        REQUIRE(species1D3V.m() == mass);
        REQUIRE(species1D3V.q() == charge);

        REQUIRE(species2D3V.m() == mass);
        REQUIRE(species2D3V.q() == charge);

        REQUIRE(species3D3V.m() == mass);
        REQUIRE(species3D3V.q() == charge);
        
        REQUIRE(species1D3V.n() == 0);
        REQUIRE(species2D3V.n() == 0);
        REQUIRE(species3D3V.n() == 0);

        species1D3V.add(5);
        species2D3V.add(5);
        species3D3V.add(5);

        REQUIRE(species1D3V.n() == 5);
        REQUIRE(species2D3V.n() == 5);
        REQUIRE(species3D3V.n() == 5);
    }

    SECTION("Add species with sampler") {
        TestSpecies1D3V species1D3V(0, 1.0);
        TestSpecies2D3V species2D3V(0, 2.0);
        TestSpecies3D3V species3D3V(0, 3.0);
        
        auto sampler1D = [](Vec<3>& v, Vec<1>& x) {
            x = {10.0};
            v = {1.0, 2.0, 3.0};
        };

        auto sampler2D = [](Vec<3>& v, Vec<2>& x) {
            x = {10.0, 20.0};
            v = {1.0, 2.0, 3.0};
        };

        auto sampler3D = [](Vec<3>& v, Vec<3>& x) {
            x = {10.0, 20.0, 30.0};
            v = {1.0, 2.0, 3.0};
        };

        species1D3V.add(1, sampler1D);
        species2D3V.add(2, sampler2D);
        species3D3V.add(3, sampler3D);


        REQUIRE(species1D3V.n() == 1);
        REQUIRE(species2D3V.n() == 2);
        REQUIRE(species3D3V.n() == 3);
        
        auto* pos_ptr1D = species1D3V.x();
        auto* vel_ptr1D = species1D3V.v();

        auto* pos_ptr2D = species2D3V.x();
        auto* vel_ptr2D = species2D3V.v();

        auto* pos_ptr3D = species3D3V.x();
        auto* vel_ptr3D = species3D3V.v();

        REQUIRE(pos_ptr1D[0].x == Catch::Approx(10.0));
        REQUIRE(vel_ptr1D[0].x == Catch::Approx(1.0));
        REQUIRE(vel_ptr1D[0].y == Catch::Approx(2.0));
        REQUIRE(vel_ptr1D[0].z == Catch::Approx(3.0));

        REQUIRE(pos_ptr2D[0].x == Catch::Approx(10.0));
        REQUIRE(pos_ptr2D[0].y == Catch::Approx(20.0));
        REQUIRE(vel_ptr2D[0].x == Catch::Approx(1.0));
        REQUIRE(vel_ptr2D[0].y == Catch::Approx(2.0));
        REQUIRE(vel_ptr2D[0].z == Catch::Approx(3.0));

        REQUIRE(pos_ptr3D[0].x == Catch::Approx(10.0));
        REQUIRE(pos_ptr3D[0].y == Catch::Approx(20.0));
        REQUIRE(pos_ptr3D[0].z == Catch::Approx(30.0));
        REQUIRE(vel_ptr3D[0].x == Catch::Approx(1.0));
        REQUIRE(vel_ptr3D[0].y == Catch::Approx(2.0));
        REQUIRE(vel_ptr3D[0].z == Catch::Approx(3.0));
    }

    SECTION("Add copy") {
        TestSpecies1D3V species1D3V(0, 1.0);
        species1D3V.add(1, [](Vec<3>& v, Vec<1>& x) {
            x = {5.5};
            v = {7.7, 8.8, 9.9};
        });

        TestSpecies2D3V species2D3V(0, 1.0);
        species2D3V.add(1, [](Vec<3>& v, Vec<2>& x) {
            x = {5.5, 6.6};
            v = {7.7, 8.8, 9.9};
        });

        TestSpecies3D3V species3D3V(0, 1.0);
        species3D3V.add(1, [](Vec<3>& v, Vec<3>& x) {
            x = {5.5, 6.6, 7.7};
            v = {7.7, 8.8, 9.9};
        });

        species1D3V.add_copy(0);
        species2D3V.add_copy(0);
        species3D3V.add_copy(0);

        REQUIRE(species1D3V.n() == 2);
        REQUIRE(species2D3V.n() == 2);
        REQUIRE(species3D3V.n() == 2);

        auto* pos_ptr1D = species1D3V.x();
        auto* vel_ptr1D = species1D3V.v();

        auto* pos_ptr2D = species2D3V.x();
        auto* vel_ptr2D = species2D3V.v();

        auto* pos_ptr3D = species3D3V.x();
        auto* vel_ptr3D = species3D3V.v();

        REQUIRE(pos_ptr1D[0] == pos_ptr1D[1]);
        REQUIRE(vel_ptr1D[0] == vel_ptr1D[1]);

        REQUIRE(pos_ptr2D[0] == pos_ptr2D[1]);
        REQUIRE(vel_ptr2D[0] == vel_ptr2D[1]);

        REQUIRE(pos_ptr3D[0] == pos_ptr3D[1]);
        REQUIRE(vel_ptr3D[0] == vel_ptr3D[1]);
    }

    SECTION("Particle removal") {
        TestSpecies1D3V species1D3V(0, 1.0);
        TestSpecies2D3V species2D3V(0, 1.0);
        TestSpecies3D3V species3D3V(0, 1.0);


        species1D3V.add(1, [](auto& v, auto& x){ x = {1.0}; });
        species1D3V.add(1, [](auto& v, auto& x){ x = {2.0}; });
        species1D3V.add(1, [](auto& v, auto& x){ x = {3.0}; });

        species2D3V.add(1, [](auto& v, auto& x){ x = {1.0, 1.0}; });
        species2D3V.add(1, [](auto& v, auto& x){ x = {2.0, 2.0}; });
        species2D3V.add(1, [](auto& v, auto& x){ x = {3.0, 3.0}; });

        species3D3V.add(1, [](auto& v, auto& x){ x = {1.0, 1.0, 1.0}; });
        species3D3V.add(1, [](auto& v, auto& x){ x = {2.0, 2.0, 2.0}; });
        species3D3V.add(1, [](auto& v, auto& x){ x = {3.0, 3.0, 3.0}; });

        REQUIRE(species1D3V.n() == 3);
        REQUIRE(species2D3V.n() == 3);
        REQUIRE(species3D3V.n() == 3);

        Vec<1> last_particle_pos1D = species1D3V.x()[2];
        Vec<2> last_particle_pos2D = species2D3V.x()[2];
        Vec<3> last_particle_pos3D = species3D3V.x()[2];
        
        species1D3V.remove(1);
        species2D3V.remove(1);
        species3D3V.remove(1);

        REQUIRE(species1D3V.n() == 2);
        REQUIRE(species2D3V.n() == 2);
        REQUIRE(species3D3V.n() == 2);

        REQUIRE(species1D3V.x()[0].x == Catch::Approx(1.0));
        REQUIRE(species2D3V.x()[0].x == Catch::Approx(1.0));
        REQUIRE(species3D3V.x()[0].x == Catch::Approx(1.0));
        
        REQUIRE(species1D3V.x()[1] == last_particle_pos1D);
        REQUIRE(species2D3V.x()[1] == last_particle_pos2D);
        REQUIRE(species3D3V.x()[1] == last_particle_pos3D);
    }
}