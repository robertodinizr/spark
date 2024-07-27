#include "kn/interpolate/field.h"

void kn::interpolate::field_at_particles(const kn::spatial::UniformGrid &field, kn::particle::ChargedSpecies &species) {
    
    size_t n = species.n();
    double* f = species.f();
    double* x = species.x();

    const double dx = field.dx();
    double* e = field.data_ptr();

    for(size_t i = 0; i < n; i++) {
        double xp = x[i];
        size_t il = static_cast<size_t>(floor(xp / dx));
        double xl = static_cast<double>(il) * dx;
        double xr = static_cast<double>(il + 1) * dx;

        f[i] = e[il]      * (xr - xp) / dx +
                e[il + 1] * (xp - xl) / dx;
    }
}