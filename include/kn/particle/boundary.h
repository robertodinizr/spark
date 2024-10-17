#pragma once

#include "kn/particle/species.h"
namespace kn::particle {

    void apply_symmetric_boundary(ChargedSpecies<1,1>& species, double xmin, double xmax);
    void apply_absorbing_boundary(ChargedSpecies<1,3>& species, double xmin, double xmax);

}