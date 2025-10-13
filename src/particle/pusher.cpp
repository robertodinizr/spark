#include "spark/particle/pusher.h"

#include "spark/particle/species.h"
#include <cmath>

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

        x[i].x += v[i].x * dt;

        const double r_old = x[i].y;

        double v_r = v[i].y;
        double v_theta = v[i].z;

        v_r = v[i].y;

        double alpha = 0.0;
        if (std::abs(r_old) > 0.0) {
            alpha = v_theta / r_old * dt;
        } else {
            v_theta = 0.0;
            v[i].z = 0.0;
            alpha = 0.0;
        }
        const double cos = std::cos(alpha);
        const double sin = std::sin(alpha);

        const double v_r_new = v_r * cos + v_theta * sin;
        const double v_theta_new = -v_r * sin + v_theta * cos;

        v[i].y = v_r_new;
        v[i].z = v_theta_new;

        x[i].y += v[i].y * dt;

    }
}

template<unsigned NX>
void spark::particle::boris_mover(spark::particle::ChargedSpecies<NX, 3>& species,
                                     const spark::core::TMatrix<core::Vec<3>, 1>& electric_field,
                                     const spark::core::TMatrix<core::Vec<3>, 1>& magnetic_field,
                                     const double dt,
                                     const bool cylindrical) {
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* E = electric_field.data_ptr();
    const auto* B = magnetic_field.data_ptr();
    const double k = (species.q() * dt) / (2.0 * species.m());

    for (size_t i = 0; i < n; i++) {
        core::Vec<3> v_minus = v[i] + E[i] * k;
        core::Vec<3> v_plus;
    
        if (B[i].norm() == 0.0) {
                v_plus = v_minus;
        } else if (B[i].norm() > 0.0) {
                const double t_mag = species.q() * B[i].norm() * dt / (2.0 * species.m());
                const core::Vec<3> t = B[i].normalized() * t_mag;
                const double s_mag = 2.0 * t_mag / (1.0 + t_mag * t_mag);
                const core::Vec<3> s = B[i].normalized() * s_mag;

                const core::Vec<3> v_prime = v_minus + cross(v_minus, t);

                v_plus = v_minus + cross(v_prime, s);
        } else {
                v_plus = v_minus;
        }

        v[i] = v_plus + E[i] * k;

        if (cylindrical) {
            if constexpr (NX >= 1) x[i].x += v[i].x * dt;

            double r_old = 0.0;
            if constexpr (NX >= 2) r_old = x[i].y;

            double v_r = v[i].y;
            double v_omega = v[i].z;

            double alpha = 0.0;
            if (r_old != 0.0) {
                alpha = v_omega / r_old * dt;
            } else {
                v_omega = 0.0;
                v[i].z = 0.0;
                alpha = 0.0;
            }

            const double cos = std::cos(alpha);
            const double sin = std::sin(alpha);

            const double v_r_new = v_r * cos + v_omega * sin;
            const double v_omega_new = -v_r * sin + v_omega * cos;

            v[i].y = v_r_new;
            v[i].z = v_omega_new;

            if constexpr (NX >= 2) x[i].y += v[i].y * dt;
        } else {
            if constexpr (NX >= 1) x[i].x += v[i].x * dt;
            if constexpr (NX >= 2) x[i].y += v[i].y * dt;
            if constexpr (NX >= 3) x[i].z += v[i].z * dt;
        }
    }
}

template void spark::particle::boris_mover<1>(spark::particle::ChargedSpecies<1, 3>&,
                                              const spark::core::TMatrix<spark::core::Vec<3>, 1>&,
                                              const spark::core::TMatrix<spark::core::Vec<3>, 1>&,
                                              const double,
                                              const bool);

template void spark::particle::boris_mover<2>(spark::particle::ChargedSpecies<2, 3>&,
                                              const spark::core::TMatrix<spark::core::Vec<3>, 1>&,
                                              const spark::core::TMatrix<spark::core::Vec<3>, 1>&,
                                              const double,
                                              const bool);
