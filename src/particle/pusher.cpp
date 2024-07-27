#include "kn/particle/pusher.h"

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