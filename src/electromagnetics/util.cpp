#include "spark/constants/constants.h"
#include "spark/electromagnetics/poisson.h"

void spark::electromagnetics::charge_density(double particle_weight,
                                             const spark::spatial::UniformGrid<1>& ion_density,
                                             const spark::spatial::UniformGrid<1>& electron_density,
                                             spark::spatial::UniformGrid<1>& out) {
    auto& out_data = out.data();
    auto& ne = electron_density.data();
    auto& ni = ion_density.data();
    double k = spark::constants::e * particle_weight / ion_density.dx().x;

    for (size_t i = 0; i < out.n().x; i++) {
        out_data(i) = k * (ni(i) - ne(i));
    }
}
