#pragma once

#include "kn/core/vec.h"
#include "kn/particle/species.h"
#include <vector>
#include <memory>
#include <unordered_set>

namespace kn::collisions {

    struct CrossSection {
        double threshold = 0.0;
        std::vector<double> energy;
        std::vector<double> cross_section;
    };

    template <unsigned NX, unsigned NV>
    class Reaction {
    public:
        CrossSection cross_section;
        virtual bool react(particle::ChargedSpecies<NX, NV>& projectile, 
                           size_t id, 
                           double kinetic_energy) = 0;
        virtual ~Reaction() {};
    };

    template <unsigned NX, unsigned NV>
    using Reactions = std::vector<std::unique_ptr<Reaction<NX, NV>>>;

    template <unsigned NX, unsigned NV> 
    class Target {
    public:
        virtual double dens_at(const core::Vec<NX>& pos) = 0;
        virtual double dens_max() = 0;
        virtual double p_null(
            const Reactions<NX, NV>& reactions,
            const particle::ChargedSpecies<NX, NV>& m_projectile) = 0; 
        virtual double temperature();
        virtual ~Target() {};
    };

    enum class RelativeDynamics {
        SlowProjectile,
        FastProjectile
    };

    template <unsigned NX, unsigned NV>
    struct ReactionConfig {
        double m_dt;
        double m_m_dx;
        std::unique_ptr<Target<NX, NV>> m_target;
        Reactions<NX, NV> m_reactions;
        RelativeDynamics m_dyn;
    };

    template <unsigned NX, unsigned NV>
    class MCCReactionSet {
    public:

        MCCReactionSet<NX, NV>(particle::ChargedSpecies<NX, NV>& projectile,
                               ReactionConfig<NX,NV>&& config);
        void react_all();

    private:

        particle::ChargedSpecies<NX, NV>& m_projectile;
        ReactionConfig<NX,NV> m_config;

        std::vector<size_t> m_particle_samples;
        std::unordered_set<size_t> m_used_cache;
        double m_max_sigmav = 0.0;  
    };



}