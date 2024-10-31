#pragma once

#include "kn/particle/species.h"

namespace kn::collisions::scattering {

double random_chi();

kn::core::Vec<3> isotropic_scatter(const kn::core::Vec<3>& v, double chi);

template <unsigned NX>
void isotropic_coll(particle::ChargedSpecies<NX, 3>& species, size_t idx, double vmag, double chi);

double electron_elastic_vmag(double kinetic_energy, double chi, double ion_mass);
double electron_excitation_vmag(double kinetic_energy, double excitation_energy);
double electron_ionization_vmag(double kinetic_energy, double ionization_energy);

}  // namespace kn::collisions::scattering