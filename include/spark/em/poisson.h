#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "spark/core/vec.h"
#include "spark/spatial/grid.h"

namespace spark::em {

class ThomasPoissonSolver1D {
public:
    ThomasPoissonSolver1D(size_t n, double dx);
    void solve(const std::vector<double>& density, std::vector<double>& out, double v0, double v1);

private:
    double m_dx;
    size_t m_n;
};

void charge_density(double particle_weight,
                    const spark::spatial::UniformGrid<1>& ion_density,
                    const spark::spatial::UniformGrid<1>& electron_density,
                    spark::spatial::UniformGrid<1>& out);

enum class CellType : uint8_t { Internal, External, BoundaryDirichlet, BoundaryNeumann };

class StructPoissonSolver2D {
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

    StructPoissonSolver2D();
    explicit StructPoissonSolver2D(const DomainProp& prop, const std::vector<Region>& regions);
    StructPoissonSolver2D(StructPoissonSolver2D&& other) noexcept;
    StructPoissonSolver2D& operator=(StructPoissonSolver2D&& other) noexcept;
    ~StructPoissonSolver2D();

    void solve(core::Matrix<2>& out, const core::Matrix<2>& rho);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class CylindricalPoissonSolver {
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

    CylindricalPoissonSolver();
    explicit CylindricalPoissonSolver(const DomainProp& prop, const std::vector<Region>& regions);
    CylindricalPoissonSolver(CylindricalPoissonSolver&& other) noexcept;
    CylindricalPoissonSolver& operator=(CylindricalPoissonSolver&& other) noexcept;
    ~CylindricalPoissonSolver();

    void solve(core::Matrix<2>& out, const core::Matrix<2>& rho);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace spark::em