#pragma once

#include "kn/particle/species.h"

namespace kn::particle {

template <unsigned NX, unsigned NV>
void move_particles(kn::particle::ChargedSpecies<NX, NV>& species, double dt);

}