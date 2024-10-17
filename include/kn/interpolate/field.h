#pragma once

#include "kn/particle/species.h"
#include "kn/spatial/grid.h"

namespace kn::interpolate {

    template <class FieldType, unsigned NX, unsigned NV>
    void field_at_particles(const FieldType& field, kn::particle::ChargedSpecies<NX, NV>& species);
}