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
                                     const spark::core::TMatrix<core::Vec<3>, 1>& electric_field,
                                     const spark::core::TMatrix<core::Vec<3>, 1>& magnetic_field,
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

void boris_mover_cylindrical(ChargedSpecies<2,3>& species,
                             const core::TMatrix<core::Vec<2>,1>& E_field,
                             double dt)
{
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const double qm = species.q() / species.m();

    for (size_t i = 0; i < n; ++i) {
	// Cartesian 
	core::Vec<3> v_cart = {v[i].z, v[i].y, v[i].x};
	core::Vec<3> E_cart = {0.0, E_field[i].y, E_field[i].x};
	core::Vec<3> B_cart = {0.0, 0.0, 0.0};

        // Mover in cartesian coordinates
        // Half-step electric acceleration
        v_cart += E_cart * (qm * dt / 2.0);

        // Full magnetic rotation
	core::Vec<3> t = B_cart * (qm * dt / 2.0);

        if (t.norm() * t.norm() > 1e-12) {
            core::Vec<3> s = t * (2.0 / (1.0 + t.norm() * t.norm()));
            core::Vec<3> v_prime = v_cart + cross(v_cart, t);
            v_cart = v_cart + cross(v_prime, s);
        } else {
	      v_cart = v_cart;  
	}

        // Second half-step electric acceleration
        v_cart += E_cart * (qm * dt / 2.0);

        // Position update in cartesian coordinates and rotate to ZR
        double z_new = x[i].x + v_cart.z * dt;
        double r_old = x[i].y;
        double x_cart_new = v_cart.z * dt;
        double y_cart_new = v_cart.y * dt; 

        double r_new = std::sqrt(x_cart_new * x_cart_new + y_cart_new * y_cart_new);
        
        double cos_alpha = 1.0;
        double sin_alpha = 0.0;

        // Rotate to ZR
	if (r_new > 1e-12) {
            double cos_alpha = x_cart_new / r_new;
            double sin_alpha = y_cart_new / r_new;
        }

        double vr_old = v_cart.y;
        double vtheta_old = v_cart.z;
        double vr_new = cos_alpha * vr_old + sin_alpha * vtheta_old;
        double vtheta_new = -sin_alpha * vr_old + cos_alpha * vtheta_old;
        
        x[i].x = z_new;
	x[i].y = r_new;
	v[i].x = v_cart.z;
	v[i].y = vr_new;
	v[i].z = vtheta_new;	
    }
}

void move_particles_cylindrical(ChargedSpecies<2, 3>& species,
                                const core::TMatrix<core::Vec<2>, 1>& E_field,
                                double dt)
{
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* f = E_field.data_ptr();
    const double k = species.q() * dt / species.m();

    for (size_t i = 0; i < n; i++) {
        v[i].x += f[i].x * k;
        v[i].y += f[i].y * k;

	double r_initial = x[i].y;
	double z_initial = x[i].x;

        x[i].x += v[i].x * dt;
        x[i].y += v[i].y * dt;

	double A = v[i].z * dt;
	double B = x[i].y;

	double r = std::sqrt((A * A) + (B * B));
	double cos = B/r;
	double sin = A/r;

	double v1 = v[i].y;
	double v2 = v[i].z;

	x[i].y = r;
	v[i].y = cos * v1 + sin * v2;
	v[i].z = -sin * v1 + cos * v2;
    }
}


template void boris_mover(ChargedSpecies<1, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);
template void boris_mover(ChargedSpecies<2, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);
template void boris_mover(ChargedSpecies<3, 3>&, const core::TMatrix<core::Vec<3>, 1>&, const core::TMatrix<core::Vec<3>, 1>&, double);

} // namespace spark::particle
