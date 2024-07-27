#pragma once

#include "kn/particle/species.h"
namespace kn::particle {

    void apply_symmetric_boundary(ChargedSpecies& species, double xmin, double xmax);

}