#include "spark/interpolate/field.h"

#include <spark/spatial/grid.h>

#include "spark/particle/species.h"

using namespace spark;
using namespace spark::core;

namespace {
    template <typename T>
    T clamp(T min, T max, T d) {
        const double t = d < min ? min : d;
        return t > max ? max : t;
    }
    }  // namespace

namespace spark::interpolate{

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
                                    core::TMatrix<T, 1>& out) {
    const size_t n = species.n();
    out.resize({n});

    auto* x = species.x();

    const auto dx = field.dx();
    const auto& f = field.data();
    const auto mdx = 1.0 / dx;
    const auto nz = field.n().x;
    const auto nr = field.n().y;

    for (size_t i = 0; i < n; ++i) {
        if (x[i].y < 1e-8) { continue; }
        const auto xp = x[i] * mdx;
        const auto xmesh_lower_left = xp.template apply<std::floor>();
        const auto idx = xmesh_lower_left.template to<size_t>();

        const double r = x[i].y;
        const double rj = idx.y * dx.y;
        const double denom = rj + 0.5 * dx.y;
        
        double f1 = 1.0;
        double f2 = 1.0;
        
        if (denom > 1e-8) {
            f1 = (rj + 0.5 * (xp.y - idx.y) * dx.y) / denom;
            f2 = (rj + 0.5 * ((xp.y - idx.y) + 1) * dx.y) / denom;
        }

        const double z_local = xp.x - idx.x;
        const double r_local = xp.y - idx.y;

        const double a1 = (1 - z_local) * (1 - r_local) * static_cast<double>(f2);
        const double a2 = z_local * (1 - r_local) * static_cast<double>(f2);
        const double a3 = z_local * r_local * static_cast<double>(f1);
        const double a4 = (1 - z_local) * r_local * static_cast<double>(f1);

        out[i] = a1 * f(idx.x, idx.y) 
                + a2 * f(idx.x + 1, idx.y) 
                + a3 * f(idx.x + 1, idx.y + 1) 
                + a4 * f(idx.x, idx.y + 1);
    }
}
}  // namespace

template <typename T, unsigned NX, unsigned NV>
void spark::interpolate::field_at_particles(const spark::spatial::TUniformGrid<T, NX>& field,
                                            const spark::particle::ChargedSpecies<NX, NV>& species,
                                            core::TMatrix<T, 1>& out) {
    spark::interpolate::field_at_particles(field, species, out);
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
template double interpolate::field_at_position(const spatial::UniformGrid<2>& field,
                                               const core::Vec<2>& pos);
template void interpolate::field_at_particles_cylindrical<double, 3>(const spatial::TUniformGrid<double, 2>& field, 
                                                                     const particle::ChargedSpecies<2, 3>& species, 
                                                                     core::TMatrix<double, 1>& out);
template void interpolate::field_at_particles_cylindrical(const spatial::TUniformGrid<Vec<2>, 2>& field,
                                                          const particle::ChargedSpecies<2, 3>& species,
                                                          TMatrix<Vec<2>, 1>& out);
