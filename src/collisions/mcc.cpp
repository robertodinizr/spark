#include "kn/collisions/mcc.h"
#include "kn/constants/constants.h"
#include "kn/particle/species.h"
#include "kn/random/random.h"

#include <sys/_types/_size_t.h>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <algorithm>

using namespace kn::collisions;

namespace {
    double interpolate_cross_section(const MonteCarloCollisions::CrossSection& cs, double energy) {
        if(energy <= cs.energy.front())
            return cs.cross_section.front();
        else if(energy >= cs.energy.back())
            return cs.cross_section.back();
        else {
            auto it = std::lower_bound(cs.energy.begin(), cs.energy.end(), energy);
            
            size_t rhs = it - cs.energy.begin();

            double x0 = cs.energy[rhs - 1];
            double x1 = cs.energy[rhs];
            double y1 = cs.cross_section[rhs - 1];
            double y0 = cs.cross_section[rhs];

            return y0 + (energy - x0) * (y1 - y0) / (x1 - x0);
        }
    }

    void sample_from_sequence(size_t n, size_t range, std::vector<size_t>& sequence, std::unordered_set<size_t>& used) {

        sequence.resize(n);
        used.clear();

        for(size_t i = 0; i < n; i++) {
            size_t num = 0;
            do {
                num = kn::random::uniform_u64() % (n + 1);
            } while(used.find(num) != used.end());
        
            used.insert(num);
            sequence[i] = num;
        }
    }

    double kinetic_energy_ev(const kn::particle::ChargedSpecies& p, size_t idx) {
        return 0.0;
    }
}

MonteCarloCollisions::MonteCarloCollisions(DomainConfig config, CrossSection&& el_cs, std::vector<CrossSection>&& exc_cs, CrossSection&& iz_cs, CrossSection&& iso_cs, CrossSection&& bs_cs) : m_config(config) {

    // TODO(lui): check how are the move assignment operators implemented by
    // default, to check if this is is being moved correctly without copy. 
    m_el_cs = std::move(el_cs);
    m_exc_cs = std::move(exc_cs);
    m_iz_cs = std::move(iz_cs);
    m_iso_cs = std::move(iso_cs);
    m_bs_cs = std::move(bs_cs);

    // Initialize the MCC parameters
    init();
}

void MonteCarloCollisions::init() {   
    m_nu_prime_e = calc_nu_prime_electrons();
    m_p_null_e = calc_p_null(m_nu_prime_e);
    m_nu_prime_i = calc_nu_prime_ions();
    m_p_null_i = calc_p_null(m_nu_prime_i);
}

double MonteCarloCollisions::calc_p_null(double nu_prime) {
    return 1.0 - std::exp(-nu_prime * m_config.m_dt);
}

double MonteCarloCollisions::total_cs_electrons(double energy) {
    double cs = 0.0;
    cs += interpolate_cross_section(m_el_cs, energy);
    cs += interpolate_cross_section(m_iz_cs, energy);
    for(const auto& exc_cs : m_exc_cs)
        cs += interpolate_cross_section(exc_cs, energy);
    return cs;
}

double MonteCarloCollisions::total_cs_ions(double energy) {
    return interpolate_cross_section(m_iso_cs, energy) + interpolate_cross_section(m_bs_cs, energy);
}

double MonteCarloCollisions::nu_prime_electrons_max(const MonteCarloCollisions::CrossSection &cs) {

        double nu_prime = 0.0;
        const double rmc = kn::constants::e / kn::constants::m_e;

        for (size_t i = 0; i < cs.energy.size(); i++) {
            double energy = cs.energy[i];
            double tcs = total_cs_electrons(energy);
            double nu = m_config.m_n_neutral * tcs * std::sqrt(2.0 * energy * rmc);
            nu_prime = std::max(nu_prime, nu);
        }

        return nu_prime;
}

double MonteCarloCollisions::nu_prime_ions_max(const MonteCarloCollisions::CrossSection &cs) {

        double nu_prime = 0.0;
        const double rmc = kn::constants::e / kn::constants::m_e;

        for (size_t i = 0; i < cs.energy.size(); i++) {
            double energy = cs.energy[i];
            double tcs = total_cs_ions(energy);
            double nu = m_config.m_n_neutral * tcs * std::sqrt(4.0 * energy * rmc);
            nu_prime = std::max(nu_prime, nu);
        }

        return nu_prime;
}

double MonteCarloCollisions::calc_nu_prime_electrons() {
    double nu_prime = 0.0;

    // TODO(lui): it's not necessary to repeat the process for all the 
    // cross sections since only the energy values are evaluated. Refactor
    // later. 
    nu_prime = std::max(nu_prime, nu_prime_electrons_max(m_el_cs));
    nu_prime = std::max(nu_prime, nu_prime_electrons_max(m_iz_cs));
    for(const auto& exc_cs : m_exc_cs)
        nu_prime = std::max(nu_prime, nu_prime_electrons_max(exc_cs));

    return nu_prime;
}

double MonteCarloCollisions::calc_nu_prime_ions() {
    // TODO(lui): same as above.
    return nu_prime_ions_max(m_iso_cs);
}

int MonteCarloCollisions::collide_electrons(particle::ChargedSpecies &electrons, particle::ChargedSpecies &ions) {

    double n_null_f = m_p_null_e * (double)electrons.n();
    size_t n_null = (size_t) std::floor(n_null_f);
    n_null = (n_null_f - (double)n_null) > random::uniform() ? n_null + 1 : n_null;

    // TODO(lui): check the performance of this sequence generation.
    sample_from_sequence(n_null, electrons.n(), m_particle_samples, m_used_cache);

    double freq_ratio_0 = 0.0;
	double freq_ratio_1 = 0.0;

    for(size_t i = 0; i < n_null; i++) {
        
        size_t p_idx = m_particle_samples[i];

        double kinetic_energy = 0.0;
        // TODO(lui): continue implementation.
    }

    return 0;
}