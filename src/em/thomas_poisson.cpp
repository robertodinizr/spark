#include "spark/constants/constants.h"
#include "spark/em/poisson.h"

using namespace spark::em;

DirichletPoissonSolver1D::DirichletPoissonSolver1D(size_t n, double dx) : m_n(n), m_dx(dx) {}

void DirichletPoissonSolver1D::solve(const std::vector<double>& density,
                                     std::vector<double>& out,
                                     double v0,
                                     double v1) {
    out.resize(m_n);
    poisson_thomas(density.data() + 1, out.data() + 1, m_n - 2, m_dx, v0, v1);
    out.front() = v0;
    out.back() = v1;
}

void DirichletPoissonSolver1D::efield(const std::vector<double>& phi, std::vector<double>& out) {
    out.resize(m_n);
    efield_extrapolate(phi.data(), out.data(), m_n, m_dx);
}

void DirichletPoissonSolver1D::poisson_thomas(const double* fin,
                                              double* yout,
                                              int n,
                                              double dx,
                                              double ylhs,
                                              double yrhs) {
    const double poisson_rhs_const = -1.0 / spark::constants::eps0;
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

void DirichletPoissonSolver1D::efield_extrapolate(const double* phi,
                                                  double* eout,
                                                  int n,
                                                  double dx) {
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