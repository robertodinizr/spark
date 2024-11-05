#pragma once

#include <sys/_types/_size_t.h>

#include <unordered_set>
#include <vector>

#include "kn/particle/species.h"
#include "reaction.h"

namespace kn::collisions {

enum class RelativeDynamics { SlowProjectile, FastProjectile };

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
    MCCReactionSet(particle::ChargedSpecies<NX, NV>& projectile,
                           ReactionConfig<NX, NV>&& config);
    void react_all();

private:
    particle::ChargedSpecies<NX, NV>& m_projectile;
    ReactionConfig<NX, NV> m_config;

    std::vector<size_t> m_particle_samples;
    std::unordered_set<size_t> m_used_cache;
    double m_max_sigma_v = 0.0;
};

}  // namespace kn::collisions