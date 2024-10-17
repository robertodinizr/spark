#include "kn/interpolate/weight.h"
#include "kn/particle/species.h"

void kn::interpolate::weight_to_grid(const kn::particle::ChargedSpecies<1,1> &species, kn::spatial::UniformGrid &out) {
    const size_t n = species.n();
    auto* x = species.x();

    out.set(0.0);

    const double dx = out.dx();
    auto& g = out.data();

    for(size_t i = 0; i < n; i++) {

        double xp = x[i].x;
        size_t il = static_cast<size_t>(floor(xp / dx));
        double xl = static_cast<double>(il) * dx;
        double xr = static_cast<double>(il + 1) * dx;

        g[il] += (xr - xp) / dx;
        g[il + 1] += (xp - xl) / dx;
    }
}

void kn::interpolate::weight_to_grid(const kn::particle::ChargedSpecies<1,3> &species, kn::spatial::UniformGrid &out) {
    const size_t n = species.n();
    auto* x = species.x();

    out.set(0.0);

    const double dx = out.dx();
    auto& g = out.data();

    for(size_t i = 0; i < n; i++) {

        double xp = x[i].x;
        size_t il = static_cast<size_t>(floor(xp / dx));
        double xl = static_cast<double>(il) * dx;
        double xr = static_cast<double>(il + 1) * dx;

        g[il] += (xr - xp) / dx;
        g[il + 1] += (xp - xl) / dx;
    }

    g.front() *= 2.0;
    g.back() *= 2.0;
}