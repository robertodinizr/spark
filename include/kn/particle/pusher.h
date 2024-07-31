#pragma once

#include "kn/particle/species.h"

namespace kn::particle {

    void move_particles(kn::particle::ChargedSpecies& species, double dt);
    void move_particles(kn::particle::ChargedSpecies1D3V& species, double dt);
}