#include "kn/interpolate/weight.h"

#include "kn/particle/species.h"
#include "kn/spatial/grid.h"

namespace {
template <unsigned NV>
void weight_to_grid(const kn::particle::ChargedSpecies<1, NV>& species,
                    kn::spatial::UniformGrid& out) {
    const size_t n = species.n();
    auto* x = species.x();

    out.set(0.0);

    const double dx = out.dx();
    auto& g = out.data();
    const double mdx = 1.0 / dx;

    for (size_t i = 0; i < n; i++) {
        const double xp_dx = x[i].x * mdx;
        const double il = floor(xp_dx);
        const size_t ils = static_cast<size_t>(il);

        g[ils] += il + 1.0 - xp_dx;
        g[ils + 1] += xp_dx - il;
    }

    g.front() *= 2.0;
    g.back() *= 2.0;
}
}  // namespace

template <class GridType, unsigned NX, unsigned NV>
void kn::interpolate::weight_to_grid(const kn::particle::ChargedSpecies<NX, NV>& species,
                                     GridType& out) {
    ::weight_to_grid(species, out);
}

template void kn::interpolate::weight_to_grid(const kn::particle::ChargedSpecies<1, 1>& species,
                                              kn::spatial::UniformGrid& out);
template void kn::interpolate::weight_to_grid(const kn::particle::ChargedSpecies<1, 3>& species,
                                              kn::spatial::UniformGrid& out);
