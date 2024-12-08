#pragma once

#include <spark/core/vec.h>

#include <memory>
#include <vector>

#include "spark/spatial/grid.h"

namespace spark::electromagnetics {

class DirichletPoissonSolver1D {
public:
    DirichletPoissonSolver1D(size_t n, double dx);
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

enum class BoundaryType { Dirichlet, Neumann, None };

class StructPoissonSolver {
public:
    struct Boundary {
        BoundaryType type_ = BoundaryType::None;
        core::Vec<2> x0, xf;
    };

    struct DomainProp {
        size_t nx = 0, ny = 0;
    };

    explicit StructPoissonSolver(const DomainProp& prop, const std::vector<Boundary>& boundaries);
    ~StructPoissonSolver();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace spark::electromagnetics