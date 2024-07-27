#include "kn/electromagnetics/poisson.h"

#include <Eigen/Sparse>
#include <memory>
#include <iostream>

using namespace kn::electromagnetics;

struct SymmetricPoissonSolver::SolverImpl {
    size_t n = 0;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int>> solver;
    Eigen::VectorXd b, x;
    Eigen::SparseMatrix<double> emat;
};

SymmetricPoissonSolver::SymmetricPoissonSolver() = default;
SymmetricPoissonSolver::~SymmetricPoissonSolver() = default;

SymmetricPoissonSolver::SymmetricPoissonSolver(SymmetricPoissonSolver&& other) noexcept {
    m_impl = std::move(other.m_impl);
}

SymmetricPoissonSolver& SymmetricPoissonSolver::operator=(SymmetricPoissonSolver &&other) noexcept {
    m_impl = std::move(other.m_impl);
    return *this;
}

SymmetricPoissonSolver::SymmetricPoissonSolver(size_t n, double dx) {
    m_impl = std::make_unique<SolverImpl>();
    m_impl->n = n;

    auto mat = Eigen::SparseMatrix<double>(n, n);

        std::vector<Eigen::Triplet<double>> triplets;
        double dx2 = dx * dx;

        m_impl->b = Eigen::VectorXd(n);

        for(size_t i = 0; i < n; i++) {
            triplets.push_back(Eigen::Triplet<double>(i, i, -2.0 / dx2));

            if(i < n - 1) {
                triplets.push_back(Eigen::Triplet<double>(i+1, i, 1.0 / dx2));
                triplets.push_back(Eigen::Triplet<double>(i, i+1, 1.0 / dx2));
            }
        }

        triplets.push_back(Eigen::Triplet<double>(n-1, 0, 1.0 / dx2));
        triplets.push_back(Eigen::Triplet<double>(0, n-1, 1.0 / dx2));

        mat.setFromTriplets(triplets.begin(), triplets.end());

        m_impl->solver.compute(mat);
        if(m_impl->solver.info() != Eigen::Success) {
            std::cout << "Solver failed to decompose matrix" <<  std::endl;
        }

        m_impl->emat = Eigen::SparseMatrix<double>(n, n);
        triplets.clear();

        for(size_t i = 0; i < n - 1; i++) {
            triplets.push_back(Eigen::Triplet<double>(i+1, i, 1.0 / (2.0 * dx)));
            triplets.push_back(Eigen::Triplet<double>(i, i+1, -1.0 / (2.0 * dx)));
        }

        triplets.push_back(Eigen::Triplet<double>(n-1, 0, 1.0 / (2.0 * dx)));
        triplets.push_back(Eigen::Triplet<double>(0, n-1, -1.0 / (2.0 * dx)));
        m_impl->emat.setFromTriplets(triplets.begin(), triplets.end());
}

void SymmetricPoissonSolver::solve(const std::vector<double> &density, std::vector<double> &out) {
    size_t s = m_impl->b.size();

    for(size_t i = 0; i < s; i++) {
        m_impl->b(i) = density[i];
    }

    m_impl->x = m_impl->solver.solve(m_impl->b);
    if(m_impl->solver.info() != Eigen::Success) {
        std::cout << "Solver failed to solve system" << std::endl;
    }

    out.assign(m_impl->x.begin(), m_impl->x.end());
}

void SymmetricPoissonSolver::grad(std::vector<double> &out) {
    auto e = (- m_impl->emat * m_impl->x).eval();
    out.assign(e.begin(), e.end());
    out[0] = out.back();
}