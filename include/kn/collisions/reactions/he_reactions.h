#pragma once

#include "kn/collisions/reaction.h"
#include "kn/collisions/scattering.h"
#include "kn/particle/species.h"
#include "kn/random/random.h"
#include "kn/constants/constants.h"

namespace kn::collisions::reactions {

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
class HeElectronIonElasticCollision : public HeCollisionBase<NX, NV> {
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
class HeElectronIonExcitationCollision : public HeCollisionBase<NX, NV> {
public:
    using HeCollisionBase<NX, NV>::HeCollisionBase;

    bool react(particle::ChargedSpecies<NX, NV>& projectile,
               size_t id,
               double kinetic_energy) override {
        if (kinetic_energy < this->m_cross_section.threshold)
            return false;

        double chi = scattering::random_chi();
        scattering::isotropic_coll(projectile, id,
                       scattering::electron_excitation_vmag(kinetic_energy, this->m_cross_section.threshold),
                       chi);
        return true;
    }
};

template <unsigned NX, unsigned NV>
class HeElectronIonIonizationCollision : public HeCollisionBase<NX, NV> {
public:
    HeElectronIonIonizationCollision(particle::ChargedSpecies<NX, NV>& ions, double t_neutral, HeCollisionConfig config, CrossSection&& cs) :
        HeCollisionBase<NX, NV>(config, std::move(cs)), ions(ions), t_neutral(t_neutral) {};

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
        double vmag = scattering::electron_ionization_vmag(kinetic_energy, this->m_cross_section.threshold);

        double chi1 = scattering::random_chi();
        scattering::isotropic_coll(projectile, id, vmag, chi1);

        // Generated electron
        double chi2 = scattering::random_chi();
        scattering::isotropic_coll(projectile, p_idx_new, vmag, chi2);

        // Generated ion
        // TODO(lui): Move from std::function to something with better performance
        ions.add(1, [event_pos, ion_mass, neutral_temperature](core::Vec<3>& v, core::Vec<NX>& x) {
            x = event_pos;
            double vtemp = std::sqrt(kn::constants::kb * neutral_temperature / ion_mass);
            v = {random::normal(0.0, vtemp), random::normal(0.0, vtemp),
                 random::normal(0.0, vtemp)};
        });

        return true;
    }

private:
    particle::ChargedSpecies<NX, NV>& ions;
    double t_neutral = 0.0;
};

}  // namespace kn::collisions::reactions