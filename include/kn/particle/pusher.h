#pragma once

#include "kn/particle/species.h"

namespace kn::particle {

    void move_particles(kn::particle::ChargedSpecies<1,1>& species, double dt);
    void move_particles(kn::particle::ChargedSpecies<1,3>& species, double dt);
}