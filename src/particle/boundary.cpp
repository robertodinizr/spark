#include "spark/particle/boundary.h"

namespace spark::particle {

template <unsigned NV>
void apply_symmetric_boundary(ChargedSpecies<1, NV>& species, double xmin, double xmax) {
    size_t n = species.n();
    auto* x = species.x();

    for (size_t i = 0; i < n; i++) {
        double& pos = x[i].x;
        if (pos < xmin) {
            double dx = xmin - pos;
            pos = xmax - dx;
        } else if (pos > xmax) {
            double dx = pos - xmax;
            pos = xmin + dx;
        }
    }
}

template <unsigned NV>
void apply_absorbing_boundary(ChargedSpecies<1, NV>& species, double xmin, double xmax) {
    auto* x = species.x();

    for (long long i = species.n() - 1; i >= 0; i--) {
        double pos = x[i].x;
        if (pos < xmin || pos > xmax) {
            species.remove(i);
        }
    }
}
}  // namespace spark::particle

template void spark::particle::apply_absorbing_boundary(ChargedSpecies<1, 1>& species,
                                                     double xmin,
                                                     double xmax);
template void spark::particle::apply_absorbing_boundary(ChargedSpecies<1, 3>& species,
                                                     double xmin,
                                                     double xmax);

template void spark::particle::apply_symmetric_boundary(ChargedSpecies<1, 1>& species,
                                                     double xmin,
                                                     double xmax);
template void spark::particle::apply_symmetric_boundary(ChargedSpecies<1, 3>& species,
                                                     double xmin,
                                                     double xmax);
