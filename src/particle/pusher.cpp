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

    for (size_t i = 0; i < n; ++i) {
        core::Vec<3> v_minus = v[i] + E[i] * k;
        core::Vec<3> v_plus;
        if (B[i].norm() > 1e-12) {
            const double tan_val = std::tan(k * B[i].norm());
            const double f = tan_val / B[i].norm();
            core::Vec<3> v_prime = v_minus + cross(v_minus, B[i]) * f;
            v_plus = v_minus + cross(v_prime, B[i]) * (2.0 * f / (1.0 + tan_val * tan_val));
        } else {
            v_plus = v_minus;
        }
        core::Vec<3> v_star = v_plus + E[i] * k;

        double r_old = x[i].y;
        double x_cart_local = r_old + v_star.y * dt;
        double y_cart_local = v_star.z * dt;

        double r_new = std::sqrt(x_cart_local * x_cart_local + y_cart_local * y_cart_local);

        double cos_alpha = 1.0;
        double sin_alpha = 0.0;
        
        if (r_new > 1e-12) {
            cos_alpha = x_cart_local / r_new;
            sin_alpha = y_cart_local / r_new;
        }

        x[i].x += v_star.x * dt;
        x[i].y = r_new;

        v[i].y =  v_star.y * cos_alpha + v_star.z * sin_alpha;
        v[i].z = -v_star.y * sin_alpha + v_star.z * cos_alpha;
        v[i].x =  v_star.x;
    }
}


template void boris_mover(ChargedSpecies<1, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);
template void boris_mover(ChargedSpecies<2, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);
template void boris_mover(ChargedSpecies<3, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);

} // namespace spark::particle
