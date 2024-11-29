#pragma once

#include "spark/particle/species.h"

namespace spark::particle {

template <unsigned NX, unsigned NV>
void move_particles(spark::particle::ChargedSpecies<NX, NV>& species, double dt);

}