#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "spark/collisions/scattering.h"
#include "spark/constants/constants.h"
#include "spark/particle/species.h"
#include "spark/random/random.h"

using namespace spark::collisions;
using namespace spark::particle;
using namespace spark::core;

TEST_CASE("Scattering Angle Tests") {
    spark::random::initialize(12345);

    SECTION("random_chi and random_chi2 return valid angles") {
        for (int i = 0; i < 100; ++i) {
            double chi = scattering::random_chi();
            REQUIRE(chi >= 0.0);
            REQUIRE(chi <= spark::constants::pi);

            double chi2 = scattering::random_chi2();
            REQUIRE(chi2 >= 0.0);
            REQUIRE(chi2 <= spark::constants::pi);
        }
    }
}

TEST_CASE("Isotropic Scattering Tests") {
    spark::random::initialize(12345);
    Vec<3> v_initial = {1.0, 2.0, 3.0};
    
    SECTION("isotropic_scatter returns a unit vector") {
        Vec<3> v_scattered = scattering::isotropic_scatter(v_initial, scattering::random_chi());
        REQUIRE(v_scattered.norm() == Catch::Approx(1.0));
    }
    
    SECTION("isotropic_coll sets the correct velocity magnitude") {
        ChargedSpecies<2, 3> species(1.0, 1.0);
        species.add(1, [&](auto& v, auto& x){ v = v_initial; });

        const double vmag = 500.0;
        scattering::isotropic_coll(species, 0, vmag, scattering::random_chi());
        
        REQUIRE(species.v()[0].norm() == Catch::Approx(vmag));
    }
}

TEST_CASE("Post-Collision Velocity Magnitude Calculations") {
    const double ke_ev = 50.0;
    const double ion_mass = 40.0 * spark::constants::amc;
    
    SECTION("Electron Elastic Collision Velocity") {
        double chi = spark::constants::pi / 2.0;
        
        // delta_E = (2 * m_e / m_ion) * (1 - cos(chi))
        // ke_final = ke_inicial * (1 - delta_E)
        // v_final = sqrt(2 * e * ke_final / m_e)
        double delta_E = (2.0 * spark::constants::m_e / ion_mass) * (1.0 - std::cos(chi));
        double ke_final_joules = ke_ev * spark::constants::e * (1.0 - delta_E);
        double expected_vmag = std::sqrt(2.0 * ke_final_joules / spark::constants::m_e);
        
        double vmag = scattering::electron_elastic_vmag(ke_ev, chi, ion_mass);
        REQUIRE(vmag == Catch::Approx(expected_vmag));
    }
    
    SECTION("Electron Excitation Collision Velocity") {
        double excitation_energy = 11.5;
        
        // ke_final = ke_inicial - E_excitation
        // v_final = sqrt(2 * e * ke_final / m_e)
        double ke_final_joules = (ke_ev - excitation_energy) * spark::constants::e;
        double expected_vmag = std::sqrt(2.0 * ke_final_joules / spark::constants::m_e);

        double vmag = scattering::electron_excitation_vmag(ke_ev, excitation_energy);
        REQUIRE(vmag == Catch::Approx(expected_vmag));
    }
    
    SECTION("Electron Ionization Collision Velocity") {
        double ionization_energy = 15.7;
        
        // ke_final_each = (ke_inicial - E_ionization) / 2
        // v_final = sqrt(2 * e * ke_final_each / m_e)
        double ke_final_joules = (ke_ev - ionization_energy) * spark::constants::e;
        double expected_vmag = std::sqrt(ke_final_joules / spark::constants::m_e);
        
        double vmag = scattering::electron_ionization_vmag(ke_ev, ionization_energy);
        REQUIRE(vmag == Catch::Approx(expected_vmag));
    }
}