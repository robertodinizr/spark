#include "kn/particle/boundary.h"

namespace kn::particle {

    void apply_symmetric_boundary(ChargedSpecies &species, double xmin, double xmax) {

        size_t n = species.n();
        double* x = species.x();

        for(size_t i = 0; i < n; i++) {
            double& pos = x[i];
            if(pos < xmin) {
                double dx = xmin - pos;
                pos = xmax - dx;            
            } else if(pos > xmax) {
                double dx = pos - xmax;
                pos = xmin + dx;
            }
        }
    }

}