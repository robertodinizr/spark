#pragma once

#include <memory>
#include <unordered_set>
#include <vector>

#include "target.h"
#include "kn/particle/species.h"

namespace kn::collisions {

struct CrossSection {
    double threshold = 0.0;
    std::vector<double> energy;
    std::vector<double> cross_section;
};

template <unsigned NX, unsigned NV>
class Reaction {
public:
    Reaction(CrossSection&& cs) : m_cross_section(cs) {}
    CrossSection m_cross_section;
    virtual bool react(particle::ChargedSpecies<NX, NV>& projectile,
                       size_t id,
                       double kinetic_energy) = 0;
    virtual ~Reaction(){};
};

template <unsigned NX, unsigned NV>
using Reactions = std::vector<std::unique_ptr<Reaction<NX, NV>>>;

}  // namespace kn::collisions