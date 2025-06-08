#pragma once

#include "spark/particle/species.h"
#include "spark/spatial/grid.h"

namespace spark::interpolate {

template <unsigned NV>
void weight_to_grid(const spark::particle::ChargedSpecies<1, NV>& species,
                    spark::spatial::UniformGrid<1>& out);
    
template <unsigned NV>
void weight_to_grid(const spark::particle::ChargedSpecies<2, NV>& species,
                    spark::spatial::UniformGrid<2>& out);

template <unsigned NV>
void weight_to_grid_cylindrical(const spark::particle::ChargedSpecies<2, NV>& species, 
                                spark::spatial::UniformGrid<2>& out);
} // namespace spark::interpolate