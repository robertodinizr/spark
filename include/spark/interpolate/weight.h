#pragma once

#include "spark/particle/species.h"

namespace spark::interpolate {

template <class GridType, unsigned NX, unsigned NV>
void weight_to_grid(const spark::particle::ChargedSpecies<NX, NV>& species, GridType& out);

template <class GridType, unsigned NX, unsigned NV>
void weight_to_grid_cylindrical(const spark::particle::ChargedSpecies<NX, NV>& species, GridType& out);

}