#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "spark/particle/pusher.h"
#include "spark/particle/species.h"
#include "spark/core/vec.h"
#include "spark/core/matrix.h"
#include "spark/constants/constants.h"

using namespace spark::particle;
using namespace spark::core;
using namespace spark::constants;

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
        x = Vec<2>{1.0, 2.0};
        v = Vec<3>{1.0, 2.0, 3.0};
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

TEST_CASE("Cylindrical Boris Mover Tests", "[particle][pusher][boris_cylindrical]") {
    
    const double mass = spark::constants::m_e;
    const double charge = -spark::constants::e; 
    const double dt = 1e-11;
    
    ChargedSpecies<2, 3> species(charge, mass);
    TMatrix<Vec<3>, 1> E_field;

    SECTION("Movimento em campo nulo") {
        species.add(1, [](auto& v, auto& x){
            x = {0.0, 1.0};
            v = {1e5, 2e5, 0.0};
        });
        E_field.resize({1}); E_field.fill({0,0,0});
        
        boris_mover_cylindrical(species, E_field, dt);

        REQUIRE(species.v()[0].x == Catch::Approx(1e5));
        REQUIRE(species.v()[0].y == Catch::Approx(2e5));
        REQUIRE(species.v()[0].z == Catch::Approx(0.0));

        REQUIRE(species.x()[0].x == Catch::Approx(1e5 * dt));
        REQUIRE(species.x()[0].y == Catch::Approx(1.0 + 2e5 * dt));
    }

    SECTION("Aceleração em campo elétrico radial (E_r)") {
        const double E_r = 1e4;

        species.add(1, [](auto& v, auto& x){
            x = Vec<2>{0.0, 1.0};
            v = Vec<3>{0.0, 0.0, 0.0};
        });
        
        E_field.resize({1}); E_field.fill(Vec<3>{0, E_r, 0});

        boris_mover_cylindrical(species, E_field, dt);
        
        double expected_vr = (charge / mass) * E_r * dt;
        REQUIRE(species.v()[0].y == Catch::Approx(expected_vr));
        REQUIRE(species.v()[0].x == Catch::Approx(0.0));
        REQUIRE(species.v()[0].z == Catch::Approx(0.0));
        REQUIRE(species.x()[0].y < 1.0);
    }
    
    SECTION("Teste de singularidade: partícula cruzando o eixo r=0") {
        species.add(1, [](auto& v, auto& x){
            x = Vec<2>{0.0, 0.1};
            v = Vec<3>{0.0, -2e6, 0.0};
        });
        E_field.resize({1}); E_field.fill(Vec<3>{0,0,0});

        double large_dt = 1e-7; 
        
        boris_mover_cylindrical(species, E_field, large_dt);

        REQUIRE(species.x()[0].y >= 0.0);
        REQUIRE(!std::isnan(species.x()[0].y));
        REQUIRE(!std::isnan(species.v()[0].y));
        REQUIRE(species.v()[0].y > 0.0);
    }
}
