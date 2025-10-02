#include "spark/interpolate/weight.h"

#include "spark/particle/species.h"
#include "spark/spatial/grid.h"

namespace {
template <unsigned NV>
void weight_to_grid(const spark::particle::ChargedSpecies<1, NV>& species,
                    spark::spatial::UniformGrid<1>& out) {
    const size_t n = species.n();
    auto* x = species.x();

    out.set(0.0);

    const double dx = out.dx().x;
    auto& g = out.data().data();
    const double mdx = 1.0 / dx;

    for (size_t i = 0; i < n; i++) {
        const double xp_dx = x[i].x * mdx;
        const double il = floor(xp_dx);
        const size_t ils = static_cast<size_t>(il);

        g[ils] += il + 1.0 - xp_dx;
        g[ils + 1] += xp_dx - il;
    }

    g.front() *= 2.0;
    g.back() *= 2.0;
}

// Template specialization for 2D - Based on Birdsall and Langdon (1991), Chapter 14 Section 14.2
template <unsigned NV>
void weight_to_grid(const spark::particle::ChargedSpecies<2, NV>& species,
                    spark::spatial::UniformGrid<2>& out) {
    const auto n = species.n();
    auto* x = species.x();

    // Cache grid is used to decrease a bit the cache misses due to
    // access of different matrix rows during the iteration.
    // TODO(lui): Remove static variable
    static auto cache_grid = spark::core::TMatrix<std::array<double, 4>, 2>();
    cache_grid.resize(out.n());
    cache_grid.fill({0, 0, 0, 0});
    out.set(0.0);

    const auto [mdx, mdy] = 1.0 / out.dx();
    const auto [nx, ny] = out.n();

    auto& grid_data = out.data();

    for (size_t i = 0; i < n; ++i) {
        const double xp = x[i].x * mdx;
        const double yp = x[i].y * mdy;

        const auto jf = floor(xp);
        const auto kf = floor(yp);

        const auto j = static_cast<size_t>(jf);
        const auto k = static_cast<size_t>(kf);

        const double x_local = xp - jf;
        const double y_local = yp - kf;

        auto& c = cache_grid(j, k);
        c[0] += (1.0 - x_local) * (1.0 - y_local);
        c[1] += x_local * (1.0 - y_local);
        c[2] += (1.0 - x_local) * y_local;
        c[3] += x_local * y_local;
    }

    for (int j = 0; j < nx - 1; j++) {
        for (int k = 0; k < ny - 1; k++) {
            auto& c = cache_grid(j, k);

            grid_data(j, k) += c[0];
            grid_data(j + 1, k) += c[1];
            grid_data(j, k + 1) += c[2];
            grid_data(j + 1, k + 1) += c[3];
        }
    }

    for (size_t j = 0; j < nx; j++) {
        grid_data(j, 0) *= 2.0;
        grid_data(j, ny - 1) *= 2.0;
    }

    for (size_t k = 0; k < ny; k++) {
        grid_data(0, k) *= 2.0;
        grid_data(nx - 1, k) *= 2.0;
    }
}
}  // namespace

template <unsigned NV>
void weight_to_grid_cylindrical(const spark::particle::ChargedSpecies<2, NV>& species,
                    spark::spatial::UniformGrid<2>& out) {
    const auto n = species.n();
    auto* x = species.x();

    static auto cache_grid = spark::core::TMatrix<std::array<double, 4>, 2>();
    cache_grid.resize(out.n());
    cache_grid.fill({0, 0, 0, 0});
    out.set(0.0);

    const auto [mdx, mdy] = 1.0 / out.dx();
    const auto [nx, ny] = out.n();

    auto& grid_data = out.data();

    for (size_t i = 0; i < n; ++i) {
        const double xp = x[i].x * mdx;
        const double yp = x[i].y * mdy;

        const auto jf = floor(xp);
        const auto kf = floor(yp);

        const auto j = static_cast<size_t>(jf);
        const auto k = static_cast<size_t>(kf);

        const double x_local = xp - jf;
        const double y_local = yp - kf;

        double rj = k * out.dx().y;

        double f1 = (rj + 0.5 * y_local * out.dx().y) / (rj + 0.5 * out.dx().y);
        double f2 = (rj + 0.5 * (y_local + 1.0) * out.dx().y) / (rj + 0.5 * out.dx().y);

        auto& c = cache_grid(j, k);
        c[0] += (1.0 - x_local) * (1.0 - y_local) * f2;
        c[1] += x_local * (1.0 - y_local) * f2;
        c[2] += (1.0 - x_local) * y_local * f1;
        c[3] += x_local * y_local * f1;
    }

    for (int j = 0; j < nx - 1; j++) {
        for (int k = 0; k < ny - 1; k++) {
            auto& c = cache_grid(j, k);

            grid_data(j, k) += c[0];
            grid_data(j + 1, k) += c[1];
            grid_data(j, k + 1) += c[2];
            grid_data(j + 1, k + 1) += c[3];
        }
    }

    for (size_t j = 0; j < nx; j++) {
        grid_data(j, 0) *= 2.0;
        grid_data(j, ny - 1) *= 2.0;
    }

    for (size_t k = 0; k < ny; k++) {
        grid_data(0, k) *= 2.0;
        grid_data(nx - 1, k) *= 2.0;
    }
}

template <class GridType, unsigned NX, unsigned NV>
void spark::interpolate::weight_to_grid(const spark::particle::ChargedSpecies<NX, NV>& species,
                                        GridType& out) {
    ::weight_to_grid(species, out);
}

template <class GridType, unsigned NX, unsigned NV>
void spark::interpolate::weight_to_grid_cylindrical(
    const spark::particle::ChargedSpecies<NX, NV>& species,
    GridType& out) {
    ::weight_to_grid_cylindrical(species, out);
}

template void spark::interpolate::weight_to_grid(
    const spark::particle::ChargedSpecies<1, 1>& species,
    spark::spatial::UniformGrid<1>& out);
template void spark::interpolate::weight_to_grid(
    const spark::particle::ChargedSpecies<1, 3>& species,
    spark::spatial::UniformGrid<1>& out);
template void spark::interpolate::weight_to_grid(
    const spark::particle::ChargedSpecies<2, 1>& species,
    spark::spatial::UniformGrid<2>& out);
template void spark::interpolate::weight_to_grid(
    const spark::particle::ChargedSpecies<2, 3>& species,
    spark::spatial::UniformGrid<2>& out);
template void spark::interpolate::weight_to_grid_cylindrical(
    const spark::particle::ChargedSpecies<2, 1>& species,
    spark::spatial::UniformGrid<2>& out);
template void spark::interpolate::weight_to_grid_cylindrical(
    const spark::particle::ChargedSpecies<2, 3>& species,
    spark::spatial::UniformGrid<2>& out);