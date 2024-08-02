#include "kn/particle/pusher.h"
#include "kn/particle/species.h"

void kn::particle::move_particles(kn::particle::ChargedSpecies &species, double dt) {
    size_t n = species.n();
    double* v = species.v();
    double* x = species.x();
    double* f = species.f();

    for(size_t i = 0; i < n; i++) {
        v[i] += f[i] * dt / 2.0;
        x[i] += v[i] * dt;
        f[i] += f[i] * dt / 2.0;
    }
}

void kn::particle::move_particles(kn::particle::ChargedSpecies1D3V &species, double dt) {
    size_t n = species.n();
    auto* v = species.v();
    double* x = species.x();
    double* f = species.f();
    double k = species.q() * dt / species.m();

    for(size_t i = 0; i < n; i++) {
        v[i].x += f[i] * k;
        x[i] += v[i].x * dt;
    }
}