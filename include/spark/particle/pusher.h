#pragma once

#include "spark/core/matrix.h"
#include "spark/core/vec.h"
#include "spark/particle/species.h"

namespace spark::particle {

template <unsigned NX, unsigned NV>
void move_particles(ChargedSpecies<NX, NV>& species,
                    const core::TMatrix<core::Vec<NX>, 1>& force,
                    double dt);

template <unsigned NX>
void boris_mover(ChargedSpecies<NX, 3>& species,
                    const core::TMatrix<core::Vec<3>, 1>& electric_field,
                    const core::TMatrix<core::Vec<3>, 1>& magnetic_field,
                    double dt);

void boris_mover_cylindrical(ChargedSpecies<2, 3>& species,
                    const core::TMatrix<core::Vec<2>,1>& E_field,
                    double dt);

void move_particles_cylindrical(ChargedSpecies<2, 3>& species,
                                const core::TMatrix<core::Vec<2>, 1>& E_field,
                                double dt);
} // namespace spark::particle
