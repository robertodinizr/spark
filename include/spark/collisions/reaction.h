#pragma once

#include <unordered_set>
#include <vector>

#include "spark/core/enum_bit_ops.h"
#include "spark/particle/species.h"

namespace spark::collisions {
struct CrossSection {
    double threshold = 0.0;
    std::vector<double> energy;
    std::vector<double> cross_section;
};

enum class ReactionOutcome : uint8_t {
    NotCollided = 0,
    Collided = 1 << 0,
    ProjectileToBeRemoved = 1 << 1
};

ENUM_CLASS_BIT_OPS(ReactionOutcome, uint8_t)

template <unsigned NX, unsigned NV>
class Reaction {
public:
    explicit Reaction(CrossSection&& cs) : m_cross_section(cs) {}
    CrossSection m_cross_section;
    virtual ReactionOutcome react(particle::ChargedSpecies<NX, NV>& projectile,
                                  size_t id,
                                  double kinetic_energy) = 0;
    virtual ~Reaction() = default;
};

template <unsigned NX, unsigned NV>
using Reactions = std::vector<std::unique_ptr<Reaction<NX, NV>>>;

}  // namespace spark::collisions