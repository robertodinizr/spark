#pragma once

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
}