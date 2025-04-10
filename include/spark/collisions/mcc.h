#pragma once

#include <unordered_set>
#include <vector>

#include "reaction.h"
#include "spark/particle/species.h"

namespace spark::collisions {

enum class RelativeDynamics { SlowProjectile, FastProjectile };

template <unsigned NX, unsigned NV>
struct ReactionConfig {
    double dt;
    double mdx;
    std::unique_ptr<Target<NX, NV>> target;
    Reactions<NX, NV> reactions;
    RelativeDynamics dyn;
};

template <unsigned NX, unsigned NV>
class MCCReactionSet {
public:
    MCCReactionSet(particle::ChargedSpecies<NX, NV>& projectile, ReactionConfig<NX, NV>&& config);
    void react_all();

private:
    particle::ChargedSpecies<NX, NV>& projectile_;
    ReactionConfig<NX, NV> config_;

    std::vector<size_t> particle_samples_;
    std::unordered_set<size_t> used_cache_;
    double max_sigma_v_ = 0.0;
};

}  // namespace spark::collisions