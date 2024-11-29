#pragma once

#include "spark/particle/species.h"

namespace spark::particle {

template <unsigned NV>
void apply_symmetric_boundary(ChargedSpecies<1, NV>& species, double xmin, double xmax);

template <unsigned NV>
void apply_absorbing_boundary(ChargedSpecies<1, NV>& species, double xmin, double xmax);

}  // namespace spark::particle