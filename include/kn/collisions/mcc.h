#pragma once

#include "kn/particle/species.h"
#include <vector>
#include <unordered_set>

namespace kn::collisions {

    class MonteCarloCollisions {
    public:
        
        struct CrossSection {
            std::vector<double> energy;
            std::vector<double> cross_section;
        };

        struct DomainConfig {
            double m_dt, m_m_dx;
            double m_t_neutral;
            double m_n_neutral;
        };

        MonteCarloCollisions(DomainConfig config, CrossSection&& el_cs, std::vector<CrossSection>&& exc_cs, CrossSection&& iz_cs, CrossSection&& iso_cs, CrossSection&& bs_cs);

        // Copy constructor and assignment operator deleted
        MonteCarloCollisions(const MonteCarloCollisions&) = delete;
        MonteCarloCollisions &operator=(const MonteCarloCollisions&) = delete;

        // Move constructor and assignment operator
        MonteCarloCollisions(MonteCarloCollisions &&other) noexcept;
        MonteCarloCollisions &operator=(MonteCarloCollisions &&other) noexcept;

        int collide_electrons(particle::ChargedSpecies& electrons, particle::ChargedSpecies& ions);

        void collide_ions(particle::ChargedSpecies& ions); 
    
    private:
        double m_nu_prime_e = 0.0, m_nu_prime_i = 0.0;
        double m_p_null_e = 0.0, m_p_null_i = 0.0;

        std::vector<size_t> m_particle_samples;
        std::unordered_set<size_t> m_used_cache;

        CrossSection m_el_cs;
        std::vector<CrossSection> m_exc_cs;
        CrossSection m_iz_cs;
        CrossSection m_iso_cs; 
        CrossSection m_bs_cs;
        
        DomainConfig m_config;

        void init();
        double calc_p_null(double nu_prime);

        double calc_nu_prime_electrons();
        double calc_nu_prime_ions();
        
        double total_cs_electrons(double energy);
        double total_cs_ions(double energy);

        double nu_prime_electrons_max(const MonteCarloCollisions::CrossSection& cs);
        double nu_prime_ions_max(const MonteCarloCollisions::CrossSection& cs);
    };

}