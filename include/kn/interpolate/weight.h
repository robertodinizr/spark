#pragma once

#include "kn/spatial/grid.h"
#include "kn/particle/species.h"

namespace kn::interpolate {

    void weight_to_grid(const kn::particle::ChargedSpecies<1,1>& species, kn::spatial::UniformGrid& out);
    void weight_to_grid(const kn::particle::ChargedSpecies<1,3>& species, kn::spatial::UniformGrid& out);
}