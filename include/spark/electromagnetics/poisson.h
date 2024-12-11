#pragma once

#include <spark/core/vec.h>

#include <functional>
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

enum class CellType { Internal, External, BoundaryDirichlet, BoundaryNeumann };

class StructPoissonSolver {
public:
    struct Region {
        CellType region_type = CellType::Internal;
        core::IntVec<2> lower_left, upper_right;
        std::function<double()> input;
    };

    struct DomainProp {
        core::IntVec<2> extents;
        core::Vec<2> dx;
    };

    explicit StructPoissonSolver(const DomainProp& prop, const std::vector<Region>& regions);
    void solve(core::Matrix<2>& out, const core::Matrix<2>& rho);
    ~StructPoissonSolver();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace spark::electromagnetics