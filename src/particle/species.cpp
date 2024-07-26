#include "kn/particle/species.h"

void kn::particle::ChargedSpecies::add(size_t np) {
    
    vs.resize(n + np);
    xs.resize(n + np);
    fs.resize(n + np);
    
    for(size_t i = n; i < n + np; i++) {
        vs[i] = 0;
        xs[i] = 0;
        fs[i] = 0;
    }

    n += np;
}

void kn::particle::ChargedSpecies::add(size_t np, std::function<void(double&, double&)> sampler) {
    
    vs.resize(n + np);
    xs.resize(n + np);
    fs.resize(n + np);
    
    for(size_t i = n; i < n + np; i++) {
        sampler(vs[i], xs[i]);
        fs[i] = 0;
    }

    n += np;
}