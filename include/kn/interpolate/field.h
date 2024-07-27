#pragma once

#include "kn/particle/species.h"
#include "kn/spatial/grid.h"

namespace kn::interpolate {

    void field_at_particles(const kn::spatial::UniformGrid& field, kn::particle::ChargedSpecies& species);

}