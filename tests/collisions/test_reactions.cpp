#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "spark/collisions/basic_reactions.h"
#include "spark/particle/species.h"
#include "spark/constants/constants.h"

using namespace spark::collisions;
using namespace spark::collisions::reactions;
using namespace spark::particle;
using namespace spark::core;

struct ReactionTestFixture {
    ChargedSpecies<2, 3> electrons;
    ChargedSpecies<2, 3> ions;
    
    double kinetic_energy_ev = 100.0;
    double initial_velocity_mag;

    ReactionTestFixture() : electrons(-spark::constants::e, spark::constants::m_e), 
                            ions(spark::constants::e, 40 * spark::constants::amc) 
    {
        double ke_joules = kinetic_energy_ev * spark::constants::e;
        initial_velocity_mag = std::sqrt(2.0 * ke_joules / electrons.m());
        
        electrons.add(1, [&](auto& v, auto& x){
            v = {initial_velocity_mag, 0.0, 0.0};
            x = {1.0, 1.0};
        });
    }

    CrossSection make_dummy_cs(double threshold) {
        return {threshold, {threshold, 200.0}, {1e-20, 1e-20}};
    }
};

TEST_CASE_METHOD(ReactionTestFixture, "ElectronElasticCollision Test", "[collisions][reactions]") {
    BasicCollisionConfig config{ions.m()};
    ElectronElasticCollision<2, 3> reaction(config, make_dummy_cs(0.0));
    
    auto outcome = reaction.react(electrons, 0, kinetic_energy_ev);
    
    REQUIRE(static_cast<bool>(outcome & ReactionOutcome::Collided));
    
    double final_vmag = electrons.v()[0].norm();
    
    // delta_E_ratio = (2 * m_e / m_ion) ~ 2.7e-5.
    // v_final = v_inicial * sqrt(1 - delta_E)
    REQUIRE(final_vmag == Catch::Approx(initial_velocity_mag));
}

TEST_CASE_METHOD(ReactionTestFixture, "ExcitationCollision Test") {
    double threshold = 15.0;
    ExcitationCollision<2, 3> reaction({}, make_dummy_cs(threshold));

    SECTION("Energy below threshold") {
        auto outcome = reaction.react(electrons, 0, 10.0);
        REQUIRE(static_cast<bool>(outcome == ReactionOutcome::NotCollided));
        REQUIRE(electrons.v()[0].norm() == Catch::Approx(initial_velocity_mag));
    }
    
    SECTION("Energy above threshold") {
        auto outcome = reaction.react(electrons, 0, kinetic_energy_ev);
        REQUIRE(static_cast<bool>(outcome & ReactionOutcome::Collided));

        double vmag_expected = scattering::electron_excitation_vmag(kinetic_energy_ev, threshold);
        REQUIRE(electrons.v()[0].norm() == Catch::Approx(vmag_expected));
    }
}

TEST_CASE_METHOD(ReactionTestFixture, "IonizationCollision Test") {
    double threshold = 20.0;
    double t_neutral = 300.0;
    IonizationCollision<2, 3> reaction(&ions, t_neutral, {}, make_dummy_cs(threshold));

    SECTION("Energy below threshold") {
        auto outcome = reaction.react(electrons, 0, 15.0);
        REQUIRE(static_cast<bool>(outcome == ReactionOutcome::NotCollided));
        REQUIRE(electrons.n() == 1);
        REQUIRE(ions.n() == 0);
    }

    SECTION("Energy above threshold") {
        auto outcome = reaction.react(electrons, 0, kinetic_energy_ev);
        REQUIRE(static_cast<bool>(outcome & ReactionOutcome::Collided));

        REQUIRE(electrons.n() == 2);
        REQUIRE(ions.n() == 1);

        double vmag_expected = scattering::electron_ionization_vmag(kinetic_energy_ev, threshold);

        REQUIRE(electrons.v()[0].norm() == Catch::Approx(vmag_expected));
        REQUIRE(electrons.v()[1].norm() == Catch::Approx(vmag_expected));
        REQUIRE(electrons.x()[1].x == electrons.x()[0].x);
        REQUIRE(ions.x()[0].x == electrons.x()[0].x);
    }
}


TEST_CASE_METHOD(ReactionTestFixture, "ChargeExchangeCollision Test") {
    ChargedSpecies<2, 3> ion_projectile(spark::constants::e, 40 * spark::constants::amc);
    ion_projectile.add(1, [&](auto& v, auto& x){ v = {1000, 0, 0}; });
    
    ChargeExchangeCollision<2, 3> reaction({}, make_dummy_cs(0.0));
    auto outcome = reaction.react(ion_projectile, 0, 10.0);
    
    REQUIRE(static_cast<bool>(outcome & ReactionOutcome::Collided));
    REQUIRE(ion_projectile.v()[0].norm() == Catch::Approx(0.0));
}