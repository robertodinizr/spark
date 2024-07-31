#pragma once

#include "kn/spatial/grid.h"
#include <vector>
#include <memory>

namespace kn::electromagnetics {

    class SymmetricPoissonSolver {
    private:
        struct SolverImpl;
        std::unique_ptr<SolverImpl> m_impl;

    public:
        SymmetricPoissonSolver();
        SymmetricPoissonSolver(size_t n, double dx);
        ~SymmetricPoissonSolver();

        // Copy constructor and assignment operator deleted
        SymmetricPoissonSolver(const SymmetricPoissonSolver &) = delete;
        SymmetricPoissonSolver &operator=(const SymmetricPoissonSolver &) = delete;

        // Move constructor and assignment operator
        SymmetricPoissonSolver(SymmetricPoissonSolver &&other) noexcept;
        SymmetricPoissonSolver &operator=(SymmetricPoissonSolver &&other) noexcept;

        void solve(const std::vector<double>& density, std::vector<double>& out);
        void grad(std::vector<double>& out);
    };

    class DirichletPoissonSolver {
    public:
        DirichletPoissonSolver(size_t n, double dx);
        void solve(const std::vector<double>& density, std::vector<double>& out, double v0, double v1);
        void efield(const std::vector<double>& phi, std::vector<double> &out);

    private:
        static void poisson_thomas(const double *fin, double *yout, int n, double dx, double ylhs, double yrhs);
        static void efield_extrapolate(const double *phi, double *eout, int n, double dx);

        double m_dx;
        size_t m_n;
    };

    void charge_density(double particle_weight, const kn::spatial::UniformGrid& ion_density, const kn::spatial::UniformGrid& electron_density, kn::spatial::UniformGrid& out);
}