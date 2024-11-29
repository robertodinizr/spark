#include "spark/interpolate/field.h"

#include "spark/particle/species.h"

namespace {
template <unsigned NV>
void field_at_particles(const spark::spatial::UniformGrid& field,
                        spark::particle::ChargedSpecies<1, NV>& species) {
    const size_t n = species.n();
    auto* f = species.f();
    auto* x = species.x();

    const double dx = field.dx();
    const double* e = field.data_ptr();
    const double mdx = 1.0 / dx;

    for (size_t i = 0; i < n; i++) {
        const double xp_dx = x[i].x * mdx;
        const double il = floor(xp_dx);
        const size_t ils = static_cast<size_t>(il);

        f[i] = {e[ils] * (il + 1.0 - xp_dx) + e[ils + 1] * (xp_dx - il)};
    }
}
}  // namespace

template <class FieldType, unsigned NX, unsigned NV>
void spark::interpolate::field_at_particles(const FieldType& field,
                                         spark::particle::ChargedSpecies<NX, NV>& species) {
    ::field_at_particles(field, species);
}

template void spark::interpolate::field_at_particles(const spark::spatial::UniformGrid& field,
                                                  spark::particle::ChargedSpecies<1, 1>& species);
template void spark::interpolate::field_at_particles(const spark::spatial::UniformGrid& field,
                                                  spark::particle::ChargedSpecies<1, 3>& species);
