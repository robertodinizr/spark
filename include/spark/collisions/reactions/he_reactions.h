#pragma once

#include "spark/collisions/reaction.h"
#include "spark/collisions/scattering.h"
#include "spark/constants/constants.h"
#include "spark/particle/species.h"
#include "spark/random/random.h"

namespace spark::collisions::reactions {

struct HeCollisionConfig {
    double he_atomic_mass = 0.0;
};

template <unsigned NX, unsigned NV>
class HeCollisionBase : public Reaction<NX, NV> {
public:
    HeCollisionConfig m_config;
    HeCollisionBase(HeCollisionConfig config, CrossSection&& cs)
        : Reaction<NX, NV>(std::move(cs)), m_config(config) {}
};

template <unsigned NX, unsigned NV>
class HeElectronElasticCollision : public HeCollisionBase<NX, NV> {
public:
    using HeCollisionBase<NX, NV>::HeCollisionBase;

    bool react(particle::ChargedSpecies<NX, NV>& projectile,
               size_t id,
               double kinetic_energy) override {
        double chi = scattering::random_chi();

        scattering::isotropic_coll(
            projectile, id,
            scattering::electron_elastic_vmag(kinetic_energy, chi, this->m_config.he_atomic_mass),
            chi);

        return true;
    }
};

template <unsigned NX, unsigned NV>
class HeExcitationCollision : public HeCollisionBase<NX, NV> {
public:
    using HeCollisionBase<NX, NV>::HeCollisionBase;

    bool react(particle::ChargedSpecies<NX, NV>& projectile,
               size_t id,
               double kinetic_energy) override {
        if (kinetic_energy < this->m_cross_section.threshold)
            return false;

        double chi = scattering::random_chi();
        scattering::isotropic_coll(
            projectile, id,
            scattering::electron_excitation_vmag(kinetic_energy, this->m_cross_section.threshold),
            chi);
        return true;
    }
};

template <unsigned NX, unsigned NV>
class HeIonizationCollision : public HeCollisionBase<NX, NV> {
public:
    HeIonizationCollision(particle::ChargedSpecies<NX, NV>& ions,
                          double t_neutral,
                          HeCollisionConfig config,
                          CrossSection&& cs)
        : HeCollisionBase<NX, NV>(config, std::move(cs)), ions(ions), t_neutral(t_neutral){};

    bool react(particle::ChargedSpecies<NX, NV>& projectile,
               size_t id,
               double kinetic_energy) override {
        if (kinetic_energy < this->m_cross_section.threshold)
            return false;

        projectile.add_copy(id);
        size_t p_idx_new = projectile.n() - 1;

        auto event_pos = projectile.x()[id];

        double ion_mass = ions.m();
        double neutral_temperature = t_neutral;
        double vmag =
            scattering::electron_ionization_vmag(kinetic_energy, this->m_cross_section.threshold);

        double chi1 = scattering::random_chi();
        scattering::isotropic_coll(projectile, id, vmag, chi1);

        // Generated electron
        double chi2 = scattering::random_chi();
        scattering::isotropic_coll(projectile, p_idx_new, vmag, chi2);

        // Generated ion
        // TODO(lui): Move from std::function to something with better performance
        ions.add(1, [event_pos, ion_mass, neutral_temperature](core::Vec<3>& v, core::Vec<NX>& x) {
            x = event_pos;
            double vtemp = std::sqrt(spark::constants::kb * neutral_temperature / ion_mass);
            v = {random::normal(0.0, vtemp), random::normal(0.0, vtemp),
                 random::normal(0.0, vtemp)};
        });

        return true;
    }

private:
    particle::ChargedSpecies<NX, NV>& ions;
    double t_neutral = 0.0;
};

template <unsigned NX, unsigned NV>
class HeIonElasticCollision : public HeCollisionBase<NX, NV> {
public:
    using HeCollisionBase<NX, NV>::HeCollisionBase;

    bool react(particle::ChargedSpecies<NX, NV>& projectile,
               size_t id,
               double kinetic_energy) override {
        double chi = scattering::random_chi2();
        double cos_chi = std::cos(chi);
        double vmag = std::sqrt(2.0 * constants::e * (kinetic_energy * cos_chi * cos_chi) /
                                projectile.m());

        scattering::isotropic_coll(projectile, id, vmag, chi);

        return true;
    }
};

template <unsigned NX, unsigned NV>
class HeIonChargeExchangeCollision : public HeCollisionBase<NX, NV> {
public:
    using HeCollisionBase<NX, NV>::HeCollisionBase;

    bool react(particle::ChargedSpecies<NX, NV>& projectile,
               size_t id,
               double kinetic_energy) override {
        projectile.v()[id] = core::Vec<3>();
        return true;
    }
};

}  // namespace spark::collisions::reactions