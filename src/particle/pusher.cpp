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

template <unsigned NX>
core::Vec<3> boris_velocity(spark::particle::ChargedSpecies<NX, 3>& species,
                                     const core::TMatrix<core::Vec<3>, 1>& electric_field,
                                     const core::TMatrix<core::Vec<3>, 1>& magnetic_field,
                                     const double dt)
{
const double k = species.q() * dt / (2.0 * species.m());

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

        if (B[i].norm() > 0.0) {
            const double f = std::tan(k * B[i].norm()) / B[i].norm();
            core::Vec<3> v_prime = v_minus + cross(v_minus, B[i]) * f;
            v_plus = v_minus + cross(v_prime, B[i]) * (2.0 * f / (1.0 + f * f * B[i].norm() * B[i].norm()));
        } else {
            v_plus = v_minus;
        }
        core::Vec<3> v_star = v_plus + E[i] * k;

        const double r = x[i].y;
        const double vr = v_star.y;
        const double vtheta = v_star.z;

        const double x_new = r + vr * dt;
        const double y_new = vtheta * dt;

        const double r_n = std::sqrt(x_new * x_new + y_new * y_new);

        double cos_alpha = 1.0;
        double sin_alpha = 0.0;

        if (r_n > 0.0) {
            cos_alpha = x_new / r_n;
            sin_alpha = y_new / r_n;
        }

        const double vr_final = vr * cos_alpha + vtheta * sin_alpha;
        const double vtheta_final = -vr * sin_alpha + vtheta * cos_alpha;

        const double vz_final = v_star.x;

        v[i] = {vz_final, vr_final, vtheta_final};

        x[i].x += vz_final * dt;
        x[i].y = r_n;  
    }

}

template void boris_mover(ChargedSpecies<1, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);
template void boris_mover(ChargedSpecies<2, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);
template void boris_mover(ChargedSpecies<3, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);

} // namespace spark::particle