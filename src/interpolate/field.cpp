#include "spark/interpolate/field.h"

#include <spark/spatial/grid.h>

#include "spark/particle/species.h"

using namespace spark;
using namespace spark::core;

namespace {

template <typename T, unsigned NV>
void field_at_particles(const spatial::TUniformGrid<T, 1>& field,
                        const particle::ChargedSpecies<1, NV>& species,
                        TMatrix<T, 1>& out) {
    const size_t n = species.n();
    out.resize({n});

    auto* x = species.x();

    const double dx = field.dx().x;
    const auto& f = field.data();
    const double mdx = 1.0 / dx;

    for (size_t i = 0; i < n; i++) {
        const double xp_dx = x[i].x * mdx;
        const double il = floor(xp_dx);
        const auto ils = static_cast<size_t>(il);

        out[i] = f[ils] * (il + 1.0 - xp_dx) + f[ils + 1] * (xp_dx - il);
    }
}

template <typename T, unsigned NV>
void field_at_particles(const spatial::TUniformGrid<T, 2>& field,
                        const particle::ChargedSpecies<2, NV>& species,
                        TMatrix<T, 1>& out) {
    const size_t n = species.n();
    out.resize({n});

    auto* x = species.x();

    const auto dx = field.dx();
    const auto& f = field.data();
    const auto mdx = 1.0 / dx;

    for (size_t i = 0; i < n; ++i) {
        const auto xp = x[i] * mdx;
        const auto xmesh_lower_left = xp.template apply<std::floor>();
        const auto xmesh_upper_right = xmesh_lower_left + 1.0;
        const auto idx = xmesh_lower_left.template to<size_t>();

        const double a1 = (xmesh_upper_right.x - xp.x) * (xmesh_upper_right.y - xp.y);
        const double a2 = (xp.x - xmesh_lower_left.x) * (xmesh_upper_right.y - xp.y);
        const double a3 = (xp.x - xmesh_lower_left.x) * (xp.y - xmesh_lower_left.y);
        const double a4 = (xmesh_upper_right.x - xp.x) * (xp.y - xmesh_lower_left.y);

        out[i] = a1 * f(idx.x, idx.y) + a2 * f(idx.x + 1, idx.y) + a3 * f(idx.x + 1, idx.y + 1) +
                 a4 * f(idx.x, idx.y + 1);
    }
}

template <typename T, unsigned NV>
void field_at_particles_cylindrical(const spatial::TUniformGrid<T, 2>& field,
                        const particle::ChargedSpecies<2, NV>& species,
                        TMatrix<T, 1>& out, 
                        const double r0) {
    const size_t n = species.n();
    out.resize({n});

    auto* x = species.x();

    const auto dx = field.dx();
    const auto& f = field.data();
    const auto mdx = 1.0 / dx;

    for (size_t i = 0; i < n; ++i) {
        const auto xp = x[i] * mdx;
        const auto xmesh_lower_left = xp.template apply<std::floor>();
        const auto xmesh_upper_right = xmesh_lower_left + 1.0;
        const auto idx = xmesh_lower_left.template to<size_t>();

        double rj = r0 + xmesh_lower_left.y * dx.y;

        double f1 = (rj + 0.5 * (xp.y - xmesh_lower_left.y) * dx.y) / (rj + 0.5 * dx.y);
        double f2 = (rj + 0.5 * (xp.y - xmesh_lower_left.y + 1) * dx.y) / (rj + 0.5 * dx.y);

        const double a1 = (xmesh_upper_right.x - xp.x) * (xmesh_upper_right.y - xp.y) * f2;
        const double a2 = (xp.x - xmesh_lower_left.x) * (xmesh_upper_right.y - xp.y) * f2;
        const double a3 = (xp.x - xmesh_lower_left.x) * (xp.y - xmesh_lower_left.y) * f1;
        const double a4 = (xmesh_upper_right.x - xp.x) * (xp.y - xmesh_lower_left.y) * f1;

        out[i] = a1 * f(idx.x, idx.y) + a2 * f(idx.x + 1, idx.y) + a3 * f(idx.x + 1, idx.y + 1) +
                 a4 * f(idx.x, idx.y + 1);
    }
}

}  // namespace

template <typename T, unsigned NX, unsigned NV>
void spark::interpolate::field_at_particles(const spark::spatial::TUniformGrid<T, NX>& field,
                                            const spark::particle::ChargedSpecies<NX, NV>& species,
                                            core::TMatrix<T, 1>& out) {
    ::field_at_particles(field, species, out);
}
template <typename T, unsigned NX, unsigned NV>
void spark::interpolate::field_at_particles_cylindrical(const spark::spatial::TUniformGrid<T, NX>& field,
                                                        const spark::particle::ChargedSpecies<NX, NV>& species,
                                                        core::TMatrix<T, 1>& out, 
                                                        const double r0) {
    ::field_at_particles_cylindrical(field, species, out, r0);
}
template <typename T, unsigned NX>
T interpolate::field_at_position(const spatial::TUniformGrid<T, NX>& field, const Vec<NX>& pos) {
    const auto& f = field.data();

    const auto xp = pos / field.dx();
    const auto xmesh_lower_left = xp.template apply<std::floor>();
    const auto xmesh_upper_right = xmesh_lower_left + 1.0;
    const auto idx = xmesh_lower_left.template to<size_t>();

    const double a1 = (xmesh_upper_right.x - xp.x) * (xmesh_upper_right.y - xp.y);
    const double a2 = (xp.x - xmesh_lower_left.x) * (xmesh_upper_right.y - xp.y);
    const double a3 = (xp.x - xmesh_lower_left.x) * (xp.y - xmesh_lower_left.y);
    const double a4 = (xmesh_upper_right.x - xp.x) * (xp.y - xmesh_lower_left.y);

    return a1 * f(idx.x, idx.y) + a2 * f(idx.x + 1, idx.y) + a3 * f(idx.x + 1, idx.y + 1) +
           a4 * f(idx.x, idx.y + 1);
}

template void interpolate::field_at_particles(const spatial::UniformGrid<1>& field,
                                              const particle::ChargedSpecies<1, 1>& species,
                                              Matrix<1>& out);
template void interpolate::field_at_particles(const spatial::UniformGrid<1>& field,
                                              const particle::ChargedSpecies<1, 3>& species,
                                              Matrix<1>& out);
template void interpolate::field_at_particles(const spatial::UniformGrid<2>& field,
                                              const particle::ChargedSpecies<2, 3>& species,
                                              Matrix<1>& out);
template void interpolate::field_at_particles(const spatial::TUniformGrid<Vec<1>, 1>& field,
                                              const particle::ChargedSpecies<1, 3>& species,
                                              TMatrix<Vec<1>, 1>& out);
template void interpolate::field_at_particles(const spatial::TUniformGrid<Vec<2>, 2>& field,
                                              const particle::ChargedSpecies<2, 3>& species,
                                              TMatrix<Vec<2>, 1>& out);

template void interpolate::field_at_particles_cylindrical(const spatial::TUniformGrid<Vec<2>, 2>& field,
                                                          const particle::ChargedSpecies<2, 3>& species,
                                                          TMatrix<Vec<2>, 1>& out, 
                                                          const double r0);

template double interpolate::field_at_position(const spatial::UniformGrid<2>& field,
                                               const core::Vec<2>& pos);
