#pragma once

#include "spark/particle/species.h"

namespace spark::interpolate {

template <class FieldType, unsigned NX, unsigned NV>
void field_at_particles(const FieldType& field, spark::particle::ChargedSpecies<NX, NV>& species);

}