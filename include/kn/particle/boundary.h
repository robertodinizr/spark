#pragma once

#include "kn/particle/species.h"
namespace kn::particle {

    template <unsigned NV>
    void apply_symmetric_boundary(ChargedSpecies<1,NV>& species, double xmin, double xmax);
    
    template <unsigned NV>
    void apply_absorbing_boundary(ChargedSpecies<1,NV>& species, double xmin, double xmax);

}