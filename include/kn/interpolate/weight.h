#pragma once

#include "kn/spatial/grid.h"
#include "kn/particle/species.h"

namespace kn::interpolate {

    void weight_to_grid(const kn::particle::ChargedSpecies& species, kn::spatial::UniformGrid& out);

}