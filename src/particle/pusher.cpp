#include "spark/particle/pusher.h"

#include "spark/particle/species.h"
#include "spark/core/vec.h"
#include <cmath>
#include <vector>
#include <array>

namespace spark::particle {

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<1, 3>& species,
                                     const core::TMatrix<core::Vec<1>, 1>& force,
                                     const double dt) {
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* f = force.data_ptr();
    const double k = species.q() * dt / species.m();

    for (size_t i = 0; i < n; i++) {
        v[i].x += f[i].x * k;
        x[i].x += v[i].x * dt;
    }
}

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<2, 3>& species,
                                     const core::TMatrix<core::Vec<2>, 1>& force,
                                     const double dt) {
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* f = force.data_ptr();
    const double k = species.q() * dt / species.m();

    for (size_t i = 0; i < n; i++) {
        v[i].x += f[i].x * k;
        v[i].y += f[i].y * k;

        x[i].x += v[i].x * dt;
        x[i].y += v[i].y * dt;
    }
}
template <unsigned NV>
void move_particles_cylindrical(
    ChargedSpecies<2, NV>& species,
    const core::TMatrix<core::Vec<2>, 1>& electric_field_rz,
    double dt)
{
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();

    const double qm = species.q() / species.m();

    for (size_t i = 0; i < n; i++) {
        const double Er = electric_field_rz[i].x;
        const double Ez = electric_field_rz[i].y;

        const double vz = v[i].x;
        const double vr = v[i].y;
        const double v_theta = v[i].z;

        const double current_r = x[i].y;

        double az, ar, atheta;

        az = qm * Ez;

        if (current_r > 0.0) {
            ar = qm * Er + v_theta * v_theta / current_r;
            atheta = -vr * v_theta / current_r;
        } else {
            ar = 0.0;
            atheta = 0.0;

            if (current_r <= 0.0) {
                x[i].y = 0.0;
                v[i].y = 0.0;
                v[i].z = 0.0;
                }
        }

        v[i].x += az * dt;
        v[i].y += ar * dt;
        v[i].z += atheta * dt;

        x[i].x += v[i].x * dt;
        x[i].y += v[i].y * dt;

        if (x[i].y <= 0.0) {
            v[i].y = 0.0;
            v[i].z = 0.0;

            x[i].x = 0.0;
            x[i].y = 0.0;
        }
    }
}

template void move_particles_cylindrical(
    ChargedSpecies<2, 3>&,
    const core::TMatrix<core::Vec<2>, 1>&,
    double dt
    );

} // namespace spark::particle