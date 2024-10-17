#pragma once

#include "kn/spatial/grid.h"
#include "kn/particle/species.h"

namespace kn::interpolate {

    template <class GridType, unsigned NX, unsigned NV>
    void weight_to_grid(const kn::particle::ChargedSpecies<NX,NV>& species, GridType& out);
}