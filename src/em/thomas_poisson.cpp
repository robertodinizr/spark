#include "spark/constants/constants.h"
#include "spark/em/poisson.h"

using namespace spark::em;

namespace {
void poisson_thomas(const double* fin, double* yout, int n, double dx, double ylhs, double yrhs) {
    constexpr double poisson_rhs_const = -1.0 / spark::constants::eps0;
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
}  // namespace

ThomasPoissonSolver1D::ThomasPoissonSolver1D(size_t n, double dx) : m_n(n), m_dx(dx) {}

void ThomasPoissonSolver1D::solve(const std::vector<double>& density,
                                  std::vector<double>& out,
                                  double v0,
                                  double v1) {
    out.resize(m_n);
    poisson_thomas(density.data() + 1, out.data() + 1, m_n - 2, m_dx, v0, v1);
    out.front() = v0;
    out.back() = v1;
}
