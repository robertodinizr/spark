#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "spark/particle/pusher.h"
#include "spark/particle/species.h"
#include "spark/core/vec.h"
#include "spark/core/matrix.h"

using namespace spark::particle;
using namespace spark::core;

using Species1D3V = ChargedSpecies<1, 3>;
using Species2D3V = ChargedSpecies<2, 3>;

template <unsigned NX>
TMatrix<Vec<NX>, 1> make_force(size_t n, const Vec<NX>& value) {
    TMatrix<Vec<NX>, 1> f({n});
    for (size_t i = 0; i < n; ++i) f(i) = value;
    return f;
}

TEST_CASE("move_particles 1D3V: position and velocity update") {
    Species1D3V species(-1.0, 2.0);
    species.add(1, [](Vec<3>& v, Vec<1>& x) {
        x = {0.0};
        v = {1.0, 2.0, 3.0};
    });
    auto force = make_force<1>(1, Vec<1>{2.0});
    double dt = 0.5;

    move_particles(species, force, dt);

    // v_new = v_old + f*q*dt/m = 1.0 + 2.0*(-1)*0.5/2 = 1.0 - 0.5 = 0.5
    // x_new = x_old + v_new*dt = 0.0 + 0.5*0.5 = 0.25
    REQUIRE(species.v()[0].x == Catch::Approx(0.5));
    REQUIRE(species.x()[0].x == Catch::Approx(0.25));
}

TEST_CASE("move_particles 2D3V: position and velocity update") {
    Species2D3V species(2.0, 4.0);
    species.add(1, [](Vec<3>& v, Vec<2>& x) {
        x = {1.0, 2.0};
        v = {1.0, 2.0, 3.0};
    });
    Vec<2> force_val{4.0, -2.0};
    auto force = make_force<2>(1, force_val);
    double dt = 0.25;

    move_particles(species, force, dt);

    // k = q*dt/m = 2*0.25/4 = 0.125
    // v.x = 1.0 + 4.0*0.125 = 1.5
    // v.y = 2.0 + (-2.0)*0.125 = 1.75
    // x.x = 1.0 + 1.5*0.25 = 1.375
    // x.y = 2.0 + 1.75*0.25 = 2.4375
    REQUIRE(species.v()[0].x == Catch::Approx(1.5));
    REQUIRE(species.v()[0].y == Catch::Approx(1.75));
    REQUIRE(species.x()[0].x == Catch::Approx(1.375));
    REQUIRE(species.x()[0].y == Catch::Approx(2.4375));
}

TEST_CASE("move_particles: zero force") {
    Species1D3V species(1.0, 1.0);
    species.add(1, [](Vec<3>& v, Vec<1>& x) {
        x = {0.0};
        v = {2.0, 0.0, 0.0};
    });
    auto force = make_force<1>(1, Vec<1>{0.0});
    double dt = 1.0;

    move_particles(species, force, dt);

    REQUIRE(species.v()[0].x == Catch::Approx(2.0));
    REQUIRE(species.x()[0].x == Catch::Approx(2.0));
}

TEST_CASE("move_particles_cylindrical: normal update") {
    Species2D3V species(1.0, 2.0);
    species.add(1, [](Vec<3>& v, Vec<2>& x) {
        x = {0.0, 2.0}; // x = z, y = r
        v = {1.0, 2.0, 3.0}; // x = vz, y = vr, z = v_theta
    });
    TMatrix<Vec<2>, 1> efield({1});
    efield(0) = {4.0, 5.0}; // Er, Ez
    double dt = 0.1;

    move_particles_cylindrical(species, efield, dt);

    // qm = 1/2 = 0.5
    // az = qm * Ez = 0.5*5 = 2.5
    // ar = qm*Er + v_theta^2/r = 0.5*4 + 9/2 = 2 + 4.5 = 6.5
    // atheta = -vr*v_theta/r = -2*3/2 = -3
    // v.x = 1.0 + 2.5*0.1 = 1.25
    // v.y = 2.0 + 6.5*0.1 = 2.65
    // v.z = 3.0 + (-3)*0.1 = 2.7
    // x.x = 0.0 + 1.25*0.1 = 0.125
    // x.y = 2.0 + 2.65*0.1 = 2.265
    REQUIRE(species.v()[0].x == Catch::Approx(1.25));
    REQUIRE(species.v()[0].y == Catch::Approx(2.65));
    REQUIRE(species.v()[0].z == Catch::Approx(2.7));
    REQUIRE(species.x()[0].x == Catch::Approx(0.125));
    REQUIRE(species.x()[0].y == Catch::Approx(2.265));
}

TEST_CASE("move_particles_cylindrical: r = 0") {
    Species2D3V species(1.0, 1.0);
    species.add(1, [](Vec<3>& v, Vec<2>& x) {
        x = {0.0, 0.0}; // r = 0
        v = {1.0, 2.0, 3.0};
    });
    TMatrix<Vec<2>, 1> efield({1});
    efield(0) = {1.0, 2.0};
    double dt = 0.1;

    move_particles_cylindrical(species, efield, dt);

    // For r == 0, ar = 0, atheta = 0, v.y = 0, v.z = 0
    REQUIRE(species.v()[0].y == Catch::Approx(0.0));
    REQUIRE(species.v()[0].z == Catch::Approx(0.0));
}

TEST_CASE("move_particles_cylindrical: r < 0") {
    Species2D3V species(1.0, 1.0);
    species.add(1, [](Vec<3>& v, Vec<2>& x) {
        x = {0.0, -0.1};
        v = {0.0, -2.0, 1.0};
    });
    TMatrix<Vec<2>, 1> efield({1});
    efield(0) = {0.0, 0.0};
    double dt = 1.0;

    //r = 0.1 + (-2.0)*1.0 = -1.9
    move_particles_cylindrical(species, efield, dt);

    REQUIRE(species.x()[0].y == Catch::Approx(0.0));
    REQUIRE(species.v()[0].y == Catch::Approx(0.0));
    REQUIRE(species.v()[0].z == Catch::Approx(0.0));
}
