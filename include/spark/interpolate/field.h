#pragma once

#include "spark/particle/species.h"
#include "spark/spatial/grid.h"

namespace spark::interpolate {

template <typename T, unsigned NV>
void field_at_particles(const spatial::TUniformGrid<T, 1>& field,
                        const particle::ChargedSpecies<1, NV>& species,
                        core::TMatrix<T, 1>& out);

template <typename T, unsigned NV>
void field_at_particles(const spatial::TUniformGrid<T, 2>& field,
                        const particle::ChargedSpecies<2, NV>& species,
                        core::TMatrix<T, 1>& out);

template <typename T, unsigned NX>
T field_at_position(const spatial::TUniformGrid<T, NX>& field, const core::Vec<NX>& pos);

template <typename T, unsigned NV>
void field_at_particles_cylindrical(const spatial::TUniformGrid<T, 2>& field,
                                    const particle::ChargedSpecies<2, NV>& species,
                                    core::TMatrix<T, 1>& out);
} // namespace spark::interpolate 