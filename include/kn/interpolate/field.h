#pragma once

#include "kn/particle/species.h"
#include "kn/spatial/grid.h"

namespace kn::interpolate {

    template <unsigned NX, unsigned NV>
    void field_at_particles(const kn::spatial::UniformGrid& field, kn::particle::ChargedSpecies<NX, NV>& species);
}