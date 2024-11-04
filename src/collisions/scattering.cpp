#include "kn/collisions/scattering.h"

#include <cmath>

#include "kn/constants/constants.h"
#include "kn/random/random.h"

using namespace kn::collisions;

double scattering::random_chi() {
    return std::acos(1.0 - 2.0 * kn::random::uniform());
}

kn::core::Vec<3> scattering::isotropic_scatter(const kn::core::Vec<3>& v, double chi) {
    const auto [x, y, z] = v.normalized();

    const double phi = 2 * kn::constants::pi * kn::random::uniform();
    const double zeta = std::acos(z);

    const double k0 = std::cos(chi);
    const double k1 = std::sin(chi) / std::sin(zeta);
    const double k2 = k1 * std::sin(phi);
    const double k3 = k1 * std::cos(phi);

    return {x * k0 + y * k2 + x * z * k3, y * k0 - x * k2 + y * z * k3,
            z * k0 - (x * x + y * y) * k3};
}

template <unsigned NX>
void scattering::isotropic_coll(particle::ChargedSpecies<NX, 3>& species,
                                size_t idx,
                                double vmag,
                                double chi) {
    auto vs = isotropic_scatter(species.v()[idx], chi);
    species.v()[idx] = {vs.x * vmag, vs.y * vmag, vs.z * vmag};
}

template void scattering::isotropic_coll<1>(particle::ChargedSpecies<1, 3>& species,
                                            size_t idx,
                                            double vmag,
                                            double chi);
template void scattering::isotropic_coll<2>(particle::ChargedSpecies<2, 3>& species,
                                            size_t idx,
                                            double vmag,
                                            double chi);
template void scattering::isotropic_coll<3>(particle::ChargedSpecies<3, 3>& species,
                                            size_t idx,
                                            double vmag,
                                            double chi);

double scattering::electron_elastic_vmag(double kinetic_energy,
                                                         double chi,
                                                         double ion_mass) {
    double delta_energy = (2.0 * kn::constants::m_e / ion_mass) * (1.0 - std::cos(chi));
    return std::sqrt(2.0 * kn::constants::e * (kinetic_energy * (1.0 - delta_energy)) /
                     kn::constants::m_e);
}

double scattering::electron_excitation_vmag(double kinetic_energy, double excitation_energy) {
    return std::sqrt(2.0 * kn::constants::e * (kinetic_energy - excitation_energy) /
                     kn::constants::m_e);
}

double scattering::electron_ionization_vmag(double kinetic_energy, double ionization_energy) {
    // No x2 because of ionization energy division
    return std::sqrt(kn::constants::e * (kinetic_energy - ionization_energy) / kn::constants::m_e);
}