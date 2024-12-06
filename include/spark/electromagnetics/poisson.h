#pragma once

#include <memory>
#include <vector>

#include "spark/spatial/grid.h"

namespace spark::electromagnetics {

class DirichletPoissonSolver {
public:
    DirichletPoissonSolver(size_t n, double dx);
    void solve(const std::vector<double>& density, std::vector<double>& out, double v0, double v1);
    void efield(const std::vector<double>& phi, std::vector<double>& out);

private:
    static void poisson_thomas(const double* fin,
                               double* yout,
                               int n,
                               double dx,
                               double ylhs,
                               double yrhs);
    static void efield_extrapolate(const double* phi, double* eout, int n, double dx);

    double m_dx;
    size_t m_n;
};

void charge_density(double particle_weight,
                    const spark::spatial::UniformGrid<1>& ion_density,
                    const spark::spatial::UniformGrid<1>& electron_density,
                    spark::spatial::UniformGrid<1>& out);
}  // namespace spark::electromagnetics