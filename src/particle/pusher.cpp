#include "spark/particle/pusher.h"

#include "spark/particle/species.h"

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<1, 1>& species, double dt) {
    size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    auto* f = species.f();

    for (size_t i = 0; i < n; i++) {
        v[i].x += f[i].x * dt / 2.0;
        x[i].x += v[i].x * dt;
        f[i].x += f[i].x * dt / 2.0;
    }
}

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<1, 3>& species, double dt) {
    size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    auto* f = species.f();
    double k = species.q() * dt / species.m();

    for (size_t i = 0; i < n; i++) {
        v[i].x += f[i].x * k;
        x[i].x += v[i].x * dt;
    }
}

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<2, 3>& species, double dt) {
    size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    auto* f = species.f();
    double k = species.q() * dt / species.m();

    for (size_t i = 0; i < n; i++) {
        v[i].x += f[i].x * k;
        v[i].y += f[i].y * k;

        x[i].x += v[i].x * dt;
        x[i].y += v[i].y * dt;
    }
}