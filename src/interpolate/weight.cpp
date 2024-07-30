#include "kn/interpolate/weight.h"
#include "kn/particle/species.h"

void kn::interpolate::weight_to_grid(const kn::particle::ChargedSpecies &species, kn::spatial::UniformGrid &out) {
    const size_t n = species.n();
    double* x = species.x();

    const double dx = out.dx();
    auto& g = out.data();

    for(size_t i = 0; i < n; i++) {

        double xp = x[i];
        size_t il = static_cast<size_t>(floor(xp / dx));
        double xl = static_cast<double>(il) * dx;
        double xr = static_cast<double>(il + 1) * dx;

        g[il] += (xr - xp) / dx;
        g[il + 1] += (xp - xl) / dx;
    }
}

void kn::interpolate::weight_to_grid(const kn::particle::ChargedSpecies1D3V &species, kn::spatial::UniformGrid &out) {
    const size_t n = species.n();
    double* x = species.x();

    const double dx = out.dx();
    auto& g = out.data();

    for(size_t i = 0; i < n; i++) {

        double xp = x[i];
        size_t il = static_cast<size_t>(floor(xp / dx));
        double xl = static_cast<double>(il) * dx;
        double xr = static_cast<double>(il + 1) * dx;

        g[il] += (xr - xp) / dx;
        g[il + 1] += (xp - xl) / dx;
    }
}