#include "spark/collisions/mcc.h"

#include <parallel_hashmap/phmap.h>

#include "spark/constants/constants.h"
#include "spark/particle/species.h"
#include "spark/random/random.h"

using namespace spark::collisions;

namespace {
auto& sample_from_sequence(const size_t n, const size_t range) {
    static phmap::flat_hash_set<size_t> cache;
    cache.clear();

    while (cache.size() < n) {
        cache.insert(spark::random::uniform(range));
    }

    return cache;
}

template <unsigned NX>
double kinetic_energy_ev(const spark::particle::ChargedSpecies<NX, 3>& p, size_t idx) {
    const auto& v = p.v()[idx];
    return 0.5 * p.m() * (v.x * v.x + v.y * v.y + v.z * v.z) / spark::constants::e;
}

double collision_frequency(const double neutral_density,
                           const double cross_section,
                           const double kinetic_energy,
                           const double mass) {
    return neutral_density * cross_section *
           std::sqrt(2.0 * spark::constants::e * kinetic_energy / mass);
}

double calc_p_null(const double nu_prime, const double dt) {
    return 1.0 - std::exp(-nu_prime * dt);
}

double interpolate_cross_section(const CrossSection& cs, const double energy) {
    if (energy <= cs.energy.front())
        return cs.cross_section.front();

    if (energy >= cs.energy.back())
        return cs.cross_section.back();

    const auto it = std::lower_bound(cs.energy.begin(), cs.energy.end(), energy);

    const size_t rhs = it - cs.energy.begin();

    const double x0 = cs.energy[rhs - 1];
    const double x1 = cs.energy[rhs];
    const double y0 = cs.cross_section[rhs - 1];
    const double y1 = cs.cross_section[rhs];

    return y0 + (energy - x0) * (y1 - y0) / (x1 - x0);
}

template <unsigned NX, unsigned NV>
double total_cs(const double energy, const Reactions<NX, NV>& reactions) {
    double cs = 0.0;
    for (const auto& reaction : reactions)
        cs += interpolate_cross_section(reaction->m_cross_section, energy);
    return cs;
}

template <unsigned NX, unsigned NV>
double max_sigmav_for_cross_section(const CrossSection& cs,
                                    const Reactions<NX, NV>& reactions,
                                    const double projectile_mass,
                                    const bool slow_projectiles) {
    double nu_prime = 0.0;
    const double rmc = (slow_projectiles ? 2.0 : 1.0) * spark::constants::e / projectile_mass;

    for (double energy : cs.energy) {
        const double tcs = total_cs<NX>(energy, reactions);
        const double nu = tcs * std::sqrt(2.0 * energy * rmc);
        nu_prime = std::max(nu_prime, nu);
    }

    return nu_prime;
}

template <unsigned NX, unsigned NV>
double max_sigmav(const Reactions<NX, NV>& reactions,
                  const double projectile_mass,
                  bool slow_projectiles) {
    double sigmav = 0.0;
    for (const auto& reaction : reactions)
        sigmav = std::max(sigmav, max_sigmav_for_cross_section(reaction->m_cross_section, reactions,
                                                               projectile_mass, slow_projectiles));
    return sigmav;
}
}  // namespace

template <unsigned NX, unsigned NV>
MCCReactionSet<NX, NV>::MCCReactionSet(particle::ChargedSpecies<NX, NV>* projectile,
                                       ReactionConfig<NX, NV>&& config)
    : projectile_(projectile), config_(std::move(config)) {
    max_sigma_v_ = max_sigmav(*config_.reactions, projectile_->m(),
                              config_.dyn == RelativeDynamics::SlowProjectile);
}

template <unsigned NX, unsigned NV>
void MCCReactionSet<NX, NV>::react_all() {
    const double nu_prime = config_.target->dens_max() * max_sigma_v_;
    const double p_null = calc_p_null(nu_prime, config_.dt);

    const double n_null_f = p_null * static_cast<double>(projectile_->n());
    auto n_null = static_cast<size_t>(std::floor(n_null_f));
    n_null = n_null_f - static_cast<double>(n_null) > random::uniform() ? n_null + 1 : n_null;

    core::Vec<3> v_random;

    const auto& samples = sample_from_sequence(n_null, projectile_->n());

    // TODO(lui): Remove this static variable and use another cache method that avoids global state.
    static std::vector<size_t> to_be_removed_cache;
    to_be_removed_cache.resize(0);
    auto& reactions = *config_.reactions;

    for (size_t p_idx : samples) {
        if (config_.dyn == RelativeDynamics::SlowProjectile) {
            const double vth =
                std::sqrt(constants::kb * config_.target->temperature() / projectile_->m());

            v_random = {random::normal() * vth, random::normal() * vth, random::normal() * vth};

            auto& vp = projectile_->v()[p_idx];
            vp.x -= v_random.x;
            vp.y -= v_random.y;
            vp.z -= v_random.z;
        }

        double kinetic_energy = kinetic_energy_ev(*projectile_, p_idx);
        const double r1 = random::uniform();

        double fr0 = 0.0;
        double fr1 = 0.0;

        double dens_n = config_.target->dens_at(projectile_->x()[p_idx]);

        for (const auto& reaction : reactions) {
            fr0 = fr1;
            fr1 += collision_frequency(
                       dens_n,
                       interpolate_cross_section(
                           reaction->m_cross_section,
                           (config_.dyn == RelativeDynamics::SlowProjectile ? 0.5 : 1.0) *
                               kinetic_energy),
                       kinetic_energy, projectile_->m()) /
                   nu_prime;

            if (r1 > fr0 && r1 <= fr1) {
                const auto outcome = reaction->react(*projectile_, p_idx, kinetic_energy);
                if (static_cast<bool>(outcome & ReactionOutcome::ProjectileToBeRemoved)) {
                    to_be_removed_cache.push_back(p_idx);
                }

                break;
            }
        }

        if (config_.dyn == RelativeDynamics::SlowProjectile) {
            auto& vp = projectile_->v()[p_idx];
            vp.x += v_random.x;
            vp.y += v_random.y;
            vp.z += v_random.z;
        }
    }

    for (const auto& remove_idx : to_be_removed_cache) {
        projectile_->remove(remove_idx);
    }
}

template class spark::collisions::MCCReactionSet<1, 3>;
template class spark::collisions::MCCReactionSet<2, 3>;
template class spark::collisions::MCCReactionSet<3, 3>;
