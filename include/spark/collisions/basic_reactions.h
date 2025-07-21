#pragma once

#include "spark/collisions/reaction.h"
#include "spark/collisions/scattering.h"
#include "spark/constants/constants.h"
#include "spark/particle/species.h"
#include "spark/random/random.h"

namespace spark::collisions::reactions {

struct BasicCollisionConfig {
    double atomic_mass = 0.0;
};

template <unsigned NX, unsigned NV>
class BasicCollision : public Reaction<NX, NV> {
public:
    BasicCollision(BasicCollisionConfig config, CrossSection&& cs)
        : Reaction<NX, NV>(std::forward<CrossSection>(cs)), m_config(config) {}

protected:
    BasicCollisionConfig m_config;
};

template <unsigned NX, unsigned NV>
class ElectronElasticCollision final : public BasicCollision<NX, NV> {
public:
    using BasicCollision<NX, NV>::BasicCollision;

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          const double kinetic_energy) override {
        double chi = scattering::random_chi();

        scattering::isotropic_coll(
            projectile, id,
            scattering::electron_elastic_vmag(kinetic_energy, chi, this->m_config.atomic_mass),
            chi);

        return ReactionOutcome::Collided;
    }
};

template <unsigned NX, unsigned NV>
class ExcitationCollision final : public BasicCollision<NX, NV> {
public:
    using BasicCollision<NX, NV>::BasicCollision;

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        if (kinetic_energy < this->m_cross_section.threshold)
            return ReactionOutcome::NotCollided;

        double chi = scattering::random_chi();
        scattering::isotropic_coll(
            projectile, id,
            scattering::electron_excitation_vmag(kinetic_energy, this->m_cross_section.threshold),
            chi);
        return ReactionOutcome::Collided;
    }
};

template <unsigned NX, unsigned NV>
class IonizationCollision final : public BasicCollision<NX, NV> {
public:
    IonizationCollision(particle::ChargedSpecies<NX, NV>* ions,
                        const double t_neutral,
                        BasicCollisionConfig config,
                        CrossSection&& cs)
        : BasicCollision<NX, NV>(config, std::forward<CrossSection>(cs)),
          ions_(ions),
          t_neutral_(t_neutral) {};

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        if (kinetic_energy < this->m_cross_section.threshold)
            return ReactionOutcome::NotCollided;

        projectile.add_copy(id);
        size_t p_idx_new = projectile.n() - 1;

        auto event_pos = projectile.x()[id];

        double ion_mass = ions_->m();
        double neutral_temperature = t_neutral_;
        double v_mag =
            scattering::electron_ionization_vmag(kinetic_energy, this->m_cross_section.threshold);

        double chi1 = scattering::random_chi();
        scattering::isotropic_coll(projectile, id, v_mag, chi1);

        // Generated electron
        double chi2 = scattering::random_chi();
        scattering::isotropic_coll(projectile, p_idx_new, v_mag, chi2);

        // Generated ion
        ions_->add(1, [event_pos, ion_mass, neutral_temperature](core::Vec<3>& v,
                                                                 core::Vec<NX>& x) {
            x = event_pos;
            const double v_th = std::sqrt(spark::constants::kb * neutral_temperature / ion_mass);
            v = {random::normal(0.0, v_th), random::normal(0.0, v_th), random::normal(0.0, v_th)};
        });

        return ReactionOutcome::Collided;
    }

private:
    spark::particle::ChargedSpecies<NX, NV>* ions_ = nullptr;
    double t_neutral_ = 0.0;
};

template <unsigned NX, unsigned NV>
class IonElasticCollision final : public BasicCollision<NX, NV> {
public:
    using BasicCollision<NX, NV>::BasicCollision;

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        const double chi = scattering::random_chi2();
        const double cos_chi = std::cos(chi);
        const double v_mag =
            std::sqrt(2.0 * constants::e * (kinetic_energy * cos_chi * cos_chi) / projectile.m());

        scattering::isotropic_coll(projectile, id, v_mag, chi);

        return ReactionOutcome::Collided;
    }
};

template <unsigned NX, unsigned NV>
class IonDiffMassesElasticCollision final : public BasicCollision<NX, NV> {
public:
    using BasicCollision<NX, NV>::BasicCollision;

    IonDiffMassesElasticCollision(double projectile_mass, double target_mass, CrossSection&& cs)
        : BasicCollision<NX, NV>(BasicCollisionConfig{}, std::forward<CrossSection>(cs)) {
        gamma_ = projectile_mass / target_mass;
    }

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        const double r = random::uniform();
        const double chi = atan(sqrt(r * (1 - r)) / ((gamma_ + 1) / 2 - r));
        const double cc = cos(chi);
        double efactor = ((cc + sqrt(cc * cc + 3)) / 3);
        efactor *= efactor;

        const double v_mag =
            std::sqrt(2.0 * constants::e * (kinetic_energy * efactor) / projectile.m());

        scattering::isotropic_coll(projectile, id, v_mag, chi);

        return ReactionOutcome::Collided;
    }

private:
    double gamma_ = 0;
};

template <unsigned NX, unsigned NV>
class ChargeExchangeCollision final : public BasicCollision<NX, NV> {
public:
    using BasicCollision<NX, NV>::BasicCollision;

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        // Zero since this is going to go back to the target ref frame
        projectile.v()[id] = core::Vec<3>();
        return ReactionOutcome::Collided;
    }
};

template <unsigned NX, unsigned NV>
class ChargeExchangeWithRemovalCollision final : public BasicCollision<NX, NV> {
public:
    ChargeExchangeWithRemovalCollision(particle::ChargedSpecies<NX, NV>* target_species,
                                       const double t_neutral,
                                       BasicCollisionConfig config,
                                       CrossSection&& cs)
        : BasicCollision<NX, NV>(config, std::forward<CrossSection>(cs)),
          target_species_(target_species),
          t_neutral_(t_neutral) {};

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        const auto event_pos = projectile.x()[id];
        const double ion_mass = target_species_->m();
        const double neutral_temperature = t_neutral_;

        target_species_->add(1, [event_pos, ion_mass, neutral_temperature](core::Vec<3>& v,
                                                                           core::Vec<NX>& x) {
            x = event_pos;
            const double v_th = std::sqrt(spark::constants::kb * neutral_temperature / ion_mass);
            v = {random::normal(0.0, v_th), random::normal(0.0, v_th), random::normal(0.0, v_th)};
        });

        return ReactionOutcome::Collided | ReactionOutcome::ProjectileToBeRemoved;
    }

private:
    particle::ChargedSpecies<NX, NV>* target_species_ = nullptr;
    double t_neutral_ = 0;
};

template <unsigned NX, unsigned NV>
class ChargeExchangeDiffCollision final : public BasicCollision<NX, NV> {
public:
    ChargeExchangeDiffCollision(particle::ChargedSpecies<NX, NV>* target_species,
                                const double t_neutral,
                                BasicCollisionConfig config,
                                CrossSection&& cs)
        : BasicCollision<NX, NV>(config, std::forward<CrossSection>(cs)),
          target_species_(target_species),
          t_neutral_(t_neutral) {};

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        auto event_pos = projectile.x()[id];
        double ion_temp = t_neutral_;
        double ion_mass = target_species_->m();

        target_species_->add(1, [event_pos, ion_mass, ion_temp](core::Vec<3>& v, core::Vec<NX>& x) {
            x = event_pos;
            const double v_th = std::sqrt(spark::constants::kb * ion_temp / ion_mass);
            v = {random::normal(0.0, v_th), random::normal(0.0, v_th), random::normal(0.0, v_th)};
        });

        return ReactionOutcome::Collided | ReactionOutcome::ProjectileToBeRemoved;
    }

private:
    spark::particle::ChargedSpecies<NX, NV>* target_species_ = nullptr;
    double t_neutral_ = 0.0;
};

template <unsigned NX, unsigned NV>
class SinkCollision final : public BasicCollision<NX, NV> {
public:
    using BasicCollision<NX, NV>::BasicCollision;

    ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        return ReactionOutcome::Collided | ReactionOutcome::ProjectileToBeRemoved;
    }
};

}  // namespace spark::collisions::reactions