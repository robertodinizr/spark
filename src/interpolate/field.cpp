#include "kn/interpolate/field.h"

#include "kn/particle/species.h"

namespace {
template <unsigned NV>
void field_at_particles(const kn::spatial::UniformGrid& field,
                        kn::particle::ChargedSpecies<1, NV>& species) {
    const size_t n = species.n();
    auto* f = species.f();
    auto* x = species.x();

    const double dx = field.dx();
    const double* e = field.data_ptr();

    for (size_t i = 0; i < n; i++) {
        const double xp = x[i].x;
        const auto il = static_cast<size_t>(floor(xp / dx));
        const double xl = static_cast<double>(il) * dx;
        const double xr = static_cast<double>(il + 1) * dx;

        f[i] = {e[il] * (xr - xp) / dx + e[il + 1] * (xp - xl) / dx};
    }
}
}  // namespace

template <class FieldType, unsigned NX, unsigned NV>
void kn::interpolate::field_at_particles(const FieldType& field,
                                         kn::particle::ChargedSpecies<NX, NV>& species) {
    ::field_at_particles(field, species);
}

template void kn::interpolate::field_at_particles(const kn::spatial::UniformGrid& field,
                                                  kn::particle::ChargedSpecies<1, 1>& species);
template void kn::interpolate::field_at_particles(const kn::spatial::UniformGrid& field,
                                                  kn::particle::ChargedSpecies<1, 3>& species);
