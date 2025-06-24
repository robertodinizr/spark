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

template<unsigned NX>
void spark::particle::boris_mover(spark::particle::ChargedSpecies<NX, 3>& species,
                                     const core::TMatrix<core::Vec<3>, 1>& electric_field,
                                     const core::TMatrix<core::Vec<3>, 1>& magnetic_field,
                                     const double dt) {
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* E = electric_field.data_ptr();
    const double k = species.q() * dt / (2.0 * species.m());

    for (size_t i = 0; i < n; i++) {
        // Half step update
        core::Vec<3> v_minus = v[i] + E[i] * k;
        // Full step rotation
        core::Vec<3> v_plus;
        const auto* B = magnetic_field.data_ptr();
        if (B[i].norm() > 0.0) {
            const double f = std::tan(k * B[i].norm()) / B[i].norm();

            core::Vec<3> v_prime = v_minus + cross(v_minus, B[i]) * f;
            v_plus = v_minus + cross(v_prime, B[i]) * (2.0 * f / (1.0 + (f * f * B[i].norm() * B[i].norm())));
        } else {
                v_plus = v_minus;
        }
        // Half step update
        v[i] = v_plus + E[i] * k;

        if constexpr (NX >= 1) x[i].x += v[i].x * dt;
        if constexpr (NX >= 2) x[i].y += v[i].y * dt;
        if constexpr (NX >= 3) x[i].z += v[i].z * dt;
    }
}

void boris_mover_cylindrical(spark::particle::ChargedSpecies<2, 3>& species,
                             const core::TMatrix<core::Vec<3>, 1>& electric_field,
                             const core::TMatrix<core::Vec<3>, 1>& magnetic_field,
                             const double dt) {
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* E = electric_field.data_ptr();
    const auto* B = magnetic_field.data_ptr();
    const double k = species.q() * dt / (2.0 * species.m());
    const double lr = 0.067;
    const double lz = 0.00157;

    for (size_t i = 0; i < n; ++i) {

        const double old_z = x[i].x;
        const double old_r = x[i].y;
        // Half step update
        core::Vec<3> v_minus = v[i] + E[i] * k;
        core::Vec<3> v_plus = v_minus;

        if (B[i].norm() > 0.0) {
            const double f = std::tan(k * B[i].norm()) / B[i].norm();
            core::Vec<3> v_prime = v_minus + cross(v_minus, B[i]) * f;
            v_plus = v_minus + cross(v_prime, B[i]) * (2.0 * f / (1.0 + f * f * B[i].norm() * B[i].norm()));
        }

        core::Vec<3> v_star = v_plus + E[i] * k;
        const double vz_star = v_star.x;
        const double vr_star = v_star.y;
        const double vtheta_star = v_star.z;

        const double r_old = x[i].y;
        const double r_new = r_old + vr_star * dt;

        const double theta_new = (r_old > 0) ? std::atan2(vtheta_star * dt, r_old) : 0.0;

        x[i].x += vz_star * dt;
        x[i].y = r_new;

        if (r_new > 1e-12) {
            const double r_ratio = r_old / r_new;
            v[i].y = vr_star * std::cos(theta_new) - vtheta_star * std::sin(theta_new);
            v[i].z = (vr_star * std::sin(theta_new) + vtheta_star * std::cos(theta_new)) * r_ratio;
        } else {
            v[i].y = vr_star;
            v[i].z = 0.0;
        }
        v[i].x = vz_star;
    }
}

template void boris_mover(ChargedSpecies<1, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);
template void boris_mover(ChargedSpecies<2, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);
template void boris_mover(ChargedSpecies<3, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);

} // namespace spark::particle
