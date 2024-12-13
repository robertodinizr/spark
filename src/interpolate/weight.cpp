#include "spark/interpolate/weight.h"

#include "spark/particle/species.h"
#include "spark/spatial/grid.h"

namespace {
template <unsigned NV>
void weight_to_grid(const spark::particle::ChargedSpecies<1, NV>& species,
                    spark::spatial::UniformGrid<1>& out) {
    const size_t n = species.n();
    auto* x = species.x();

    out.set(0.0);

    const double dx = out.dx().x;
    auto& g = out.data().data();
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
void spark::interpolate::weight_to_grid(const spark::particle::ChargedSpecies<NX, NV>& species,
                                        GridType& out) {
    ::weight_to_grid(species, out);
}

template void spark::interpolate::weight_to_grid(
    const spark::particle::ChargedSpecies<1, 1>& species,
    spark::spatial::UniformGrid<1>& out);
template void spark::interpolate::weight_to_grid(
    const spark::particle::ChargedSpecies<1, 3>& species,
    spark::spatial::UniformGrid<1>& out);
