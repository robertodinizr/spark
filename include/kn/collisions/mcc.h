#pragma once

#include "kn/particle/species.h"
#include <sys/_types/_size_t.h>
#include <vector>
#include <unordered_set>

namespace kn::collisions {

    class MonteCarloCollisions {
    public:

        enum class CollisionProjectile {
            Ion,
            Electron
        };

        enum class CollisionType { 
            Ionization,
            Excitation,
            Elastic,
            Isotropic,
            Backscattering
        };
        
        struct CollisionReaction {
            double energy_threshold = 0.0;
            std::vector<double> energy;
            std::vector<double> cross_section;
            CollisionType type;
            CollisionProjectile projectile;
        };

        struct DomainConfig {
            double m_dt, m_m_dx;
            double m_t_neutral;
            double m_n_neutral;
            double m_m_ion;
        };

        MonteCarloCollisions(DomainConfig config, std::vector<CollisionReaction>&& cs);

        // Copy constructor and assignment operator deleted
        MonteCarloCollisions(const MonteCarloCollisions&) = delete;
        MonteCarloCollisions &operator=(const MonteCarloCollisions&) = delete;

        // Move constructor and assignment operator
        MonteCarloCollisions(MonteCarloCollisions &&other) noexcept;
        MonteCarloCollisions &operator=(MonteCarloCollisions &&other) noexcept;

        int collide_electrons(particle::ChargedSpecies<1, 3> &electrons,
                              particle::ChargedSpecies<1, 3> &ions);

        void collide_ions(particle::ChargedSpecies<1, 3> &ions);

      private:
        double m_nu_prime_e = 0.0, m_nu_prime_i = 0.0;
        double m_p_null_e = 0.0, m_p_null_i = 0.0;

        std::vector<size_t> m_particle_samples;
        std::unordered_set<size_t> m_used_cache;

        std::vector<CollisionReaction> m_electron_cs;
        std::vector<CollisionReaction> m_ion_cs;
        
        DomainConfig m_config;

        void init();
        double calc_p_null(double nu_prime);

        double calc_nu_prime_electrons();
        double calc_nu_prime_ions();
        
        double total_cs_electrons(double energy);
        double total_cs_ions(double energy);

        double nu_prime_electrons_max(const MonteCarloCollisions::CollisionReaction& cs);
        double nu_prime_ions_max(const MonteCarloCollisions::CollisionReaction& cs);

        double frequency_ratio(const MonteCarloCollisions::CollisionReaction& cs, double kinetic_energy);

        void isotropic_coll(particle::ChargedSpecies<1, 3>& species, size_t idx, double vmag, double chi);

        bool electron_elastic_coll(particle::ChargedSpecies<1, 3> &electrons,
                                     particle::ChargedSpecies<1, 3> &ions,
                                     size_t p_idx, double kinetic_energy);
        
        bool electron_excitation_coll(particle::ChargedSpecies<1, 3> &electrons,
                                     particle::ChargedSpecies<1, 3> &ions,
                                     size_t p_idx, double kinetic_energy, double threshold);

        bool electron_ionization_coll(particle::ChargedSpecies<1, 3> &electrons,
                                     particle::ChargedSpecies<1, 3> &ions,
                                     size_t p_idx, double kinetic_energy, double threshold);

        bool ions_isotropic_coll(particle::ChargedSpecies<1, 3> &ions, size_t p_idx,
                       double kinetic_energy_rel);
    };

}