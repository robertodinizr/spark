#include "spark/particle/pusher.h"

#include "spark/particle/species.h"

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

template <>
void spark::particle::move_particles_cylindrical(spark::particle::ChargedSpecies<2, 3>& species,
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

        // Axial update
        x[i].x += v[i].x * dt;

        // Save old radius to compute angular increment. Positions use (axial, radius)
        const double r_old = x[i].y;

        // Update radial velocity with radial force
        double vr = v[i].y;   // radial velocity
        double vphi = v[i].z; // azimuthal (tangential) velocity

        // Apply force to radial velocity (already done above), use updated vr
        vr = v[i].y; // v[i].y already incremented by f*y * k

        // Compute angular displacement dtheta = vphi / r * dt
        // Guard against r == 0 -> set dtheta = 0 and clear vphi to avoid singularity
        const double eps = 1e-12;
        double dtheta = 0.0;
        if (std::abs(r_old) > eps) {
            dtheta = vphi / r_old * dt;
        } else {
            // Particle at axis: angular motion undefined. Zero the azimuthal component to keep stable.
            vphi = 0.0;
            v[i].z = 0.0;
            dtheta = 0.0;
        }

        // Rotate (vr, vphi) by dtheta to account for basis rotation in cylindrical coords
        const double c = std::cos(dtheta);
        const double s = std::sin(dtheta);

        const double vr_new = vr * c - vphi * s;
        const double vphi_new = vr * s + vphi * c;

        v[i].y = vr_new;
        v[i].z = vphi_new;

        // Update radius using the rotated radial velocity
        x[i].y += v[i].y * dt;

    }
}