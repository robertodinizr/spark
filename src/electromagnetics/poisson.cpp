#include "kn/electromagnetics/poisson.h"

#include <Eigen/Sparse>
#include <iostream>
#include <memory>

#include "kn/constants/constants.h"

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

SymmetricPoissonSolver& SymmetricPoissonSolver::operator=(SymmetricPoissonSolver&& other) noexcept {
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

    for (size_t i = 0; i < n; i++) {
        triplets.push_back(Eigen::Triplet<double>(i, i, -2.0 / dx2));

        if (i < n - 1) {
            triplets.push_back(Eigen::Triplet<double>(i + 1, i, 1.0 / dx2));
            triplets.push_back(Eigen::Triplet<double>(i, i + 1, 1.0 / dx2));
        }
    }

    triplets.push_back(Eigen::Triplet<double>(n - 1, 0, 1.0 / dx2));
    triplets.push_back(Eigen::Triplet<double>(0, n - 1, 1.0 / dx2));

    mat.setFromTriplets(triplets.begin(), triplets.end());

    m_impl->solver.compute(mat);
    if (m_impl->solver.info() != Eigen::Success) {
        std::cout << "Solver failed to decompose matrix" << std::endl;
    }

    m_impl->emat = Eigen::SparseMatrix<double>(n, n);
    triplets.clear();

    for (size_t i = 0; i < n - 1; i++) {
        triplets.push_back(Eigen::Triplet<double>(i + 1, i, 1.0 / (2.0 * dx)));
        triplets.push_back(Eigen::Triplet<double>(i, i + 1, -1.0 / (2.0 * dx)));
    }

    triplets.push_back(Eigen::Triplet<double>(n - 1, 0, 1.0 / (2.0 * dx)));
    triplets.push_back(Eigen::Triplet<double>(0, n - 1, -1.0 / (2.0 * dx)));
    m_impl->emat.setFromTriplets(triplets.begin(), triplets.end());
}

void SymmetricPoissonSolver::solve(const std::vector<double>& density, std::vector<double>& out) {
    size_t s = m_impl->b.size();

    for (size_t i = 0; i < s; i++) {
        m_impl->b(i) = density[i];
    }

    m_impl->x = m_impl->solver.solve(m_impl->b);
    if (m_impl->solver.info() != Eigen::Success) {
        std::cout << "Solver failed to solve system" << std::endl;
    }

    out.assign(m_impl->x.begin(), m_impl->x.end());
}

void SymmetricPoissonSolver::grad(std::vector<double>& out) {
    auto e = (-m_impl->emat * m_impl->x).eval();
    out.assign(e.begin(), e.end());
    out[0] = out.back();
}

DirichletPoissonSolver::DirichletPoissonSolver(size_t n, double dx) : m_n(n), m_dx(dx) {}

void DirichletPoissonSolver::solve(const std::vector<double>& density,
                                   std::vector<double>& out,
                                   double v0,
                                   double v1) {
    out.resize(m_n);
    poisson_thomas(density.data() + 1, out.data() + 1, m_n - 2, m_dx, v0, v1);
    out.front() = v0;
    out.back() = v1;
}

void DirichletPoissonSolver::efield(const std::vector<double>& phi, std::vector<double>& out) {
    out.resize(m_n);
    efield_extrapolate(phi.data(), out.data(), m_n, m_dx);
}

void DirichletPoissonSolver::poisson_thomas(const double* fin,
                                            double* yout,
                                            int n,
                                            double dx,
                                            double ylhs,
                                            double yrhs) {
    const double poisson_rhs_const = -1.0 / kn::constants::eps0;
    const double k = poisson_rhs_const * dx * dx;
    double cprime = -0.5;

    yout[0] = (fin[0] * k - ylhs) / -2.0;

    for (int i = 1; i < n - 1; ++i) {
        yout[i] = (fin[i] * k - yout[i - 1]) / (-2.0 - cprime);
        cprime = 1.0 / (-2.0 - cprime);
    }

    yout[n - 1] = ((fin[n - 1] * k - yrhs) - yout[n - 2]) / (-2.0 - cprime);

    for (int i = n - 2; i >= 1; --i) {
        yout[i] -= cprime * yout[i + 1];
        cprime = -2.0 - 1.0 / cprime;
    }

    yout[0] -= -0.5 * yout[1];
}

void DirichletPoissonSolver::efield_extrapolate(const double* phi, double* eout, int n, double dx) {
    const double k = -1.0 / (2.0 * dx);

    for (int i = 1; i < n - 1; ++i) {
        eout[i] = k * (phi[i + 1] - phi[i - 1]);
    }

    // TODO(lui): check if boundaries are correctly implemented for the benchmark.
    eout[0] = 2.0 * eout[1] - eout[2];
    eout[n - 1] = 2.0 * eout[n - 2] - eout[n - 3];
    // eout[0]   = - (phi[1]   - phi[0])   / dx;
    // eout[n-1] = - (phi[n-1] - phi[n-2]) / dx;
}

void kn::electromagnetics::charge_density(double particle_weight,
                                          const kn::spatial::UniformGrid& ion_density,
                                          const kn::spatial::UniformGrid& electron_density,
                                          kn::spatial::UniformGrid& out) {
    auto& out_data = out.data();
    auto& ne = electron_density.data();
    auto& ni = ion_density.data();
    double k = kn::constants::e * particle_weight / ion_density.dx();

    for (size_t i = 0; i < out.n(); i++) {
        out_data[i] = k * (ni[i] - ne[i]);
    }
}