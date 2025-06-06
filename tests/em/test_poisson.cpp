#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "spark/em/poisson.h"
#include "spark/core/matrix.h"
#include "spark/core/vec.h"
#include <vector>
#include <cmath>

using namespace spark::em;
using namespace spark::core;

TEST_CASE("ThomasPoissonSolver1D: Initialization and zero density") {

    // phi(x) = Ax + B
    // phi(0) = B = v0
    // phi(L) = v1 = AL + v0
    // A = (v1 - v0) / L
    // phi(x) = (v1 - v0) x / L + v0

    size_t n = 5;
    double dx = 1.0;
    ThomasPoissonSolver1D solver(n, dx);

    std::vector<double> density(n, 0.0);
    std::vector<double> potential;
    double v0 = 1.0, v1 = 2.0;

    solver.solve(density, potential, v0, v1);

    for (size_t i = 0; i < n; ++i) {
        double expected = v0 + (v1 - v0) * (double(i) / (n - 1));
        REQUIRE(potential[i] == Catch::Approx(expected));
    }
}

TEST_CASE("ThomasPoissonSolver1D: Constant density") {

    // d2 phi / dx2 = - rho / e0
    // d phi / dx = - rho * x / e0 + C1
    // phi (x) = - rho * x2 / 2 * e0 + C1 * x + C2
    // phi (0) = C2 = v0
    // phi (L) = - rho * L2 / 2 * e0 + C1 * L + C2   
    // C1 = ((v1 - v0) / L) + rho * L / e0 
    // phi (x) = (- rho * x^2 / 2e0) + ((v1 - v0) * x / L) + (rho * L * x / e0) + v0

    size_t n = 5;
    double dx = 1.0;
    ThomasPoissonSolver1D solver(n, dx);

    std::vector<double> density(n, 1.0);
    std::vector<double> potential;
    double v0 = 0.0, v1 = 0.0;

    solver.solve(density, potential, v0, v1);

    double L = (n - 1) * dx;
    double factor = 1.0 / spark::constants::eps0 / 2.0;
    for (size_t i = 0; i < n; ++i) {
        double x = i * dx;
        double expected = factor * x * (L - x);
        REQUIRE(potential[i] == Catch::Approx(expected).margin(1e-10));
    }
}

TEST_CASE("StructPoissonSolver2D: Initialization and zero density, Dirichlet boundaries") {
    StructPoissonSolver2D::DomainProp prop;
    prop.extents = IntVec<2>{3, 3};
    prop.dx = Vec<2>{1.0, 1.0};

    std::vector<StructPoissonSolver2D::Region> regions = {
        {CellType::BoundaryDirichlet, {0, 0}, {0, 2}, []{ return 1.0; }},
        {CellType::BoundaryDirichlet, {2, 0}, {2, 2}, []{ return 2.0; }},
        {CellType::BoundaryDirichlet, {0, 0}, {2, 0}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {0, 2}, {2, 2}, []{ return 3.0; }},
    };

    StructPoissonSolver2D solver(prop, regions);

    Matrix<2> rho({3});
    rho.fill(0.0);
    Matrix<2> phi;

    solver.solve(phi, rho);

    REQUIRE(phi(0, 1) == Catch::Approx(1.0)); 
    REQUIRE(phi(2, 1) == Catch::Approx(2.0));
    REQUIRE(phi(1, 0) == Catch::Approx(0.0));
    REQUIRE(phi(1, 2) == Catch::Approx(3.0));

    REQUIRE(phi(0, 0) == Catch::Approx(0.0));
    REQUIRE(phi(2, 0) == Catch::Approx(0.0));
    REQUIRE(phi(0, 2) == Catch::Approx(3.0));
    REQUIRE(phi(2, 2) == Catch::Approx(3.0));
}

TEST_CASE("StructPoissonSolver2D: Constant density, Dirichlet boundaries") {
    StructPoissonSolver2D::DomainProp prop;
    prop.extents = IntVec<2>{3, 3};
    prop.dx = Vec<2>{1.0, 1.0};

    std::vector<StructPoissonSolver2D::Region> regions = {
        {CellType::BoundaryDirichlet, {0, 0}, {0, 2}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {2, 0}, {2, 2}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {0, 0}, {2, 0}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {0, 2}, {2, 2}, []{ return 0.0; }},
    };

    StructPoissonSolver2D solver(prop, regions);

    Matrix<2> rho({3});
    rho.fill(1.0);
    Matrix<2> phi;

    solver.solve(phi, rho);

    double center = phi(1, 1);
    for (int i = 1; i < 2; ++i)
        for (int j = 1; j < 2; ++j)
            REQUIRE(phi(i, j) == Catch::Approx(center));
}

TEST_CASE("CylindricalPoissonSolver2D: Initialization and zero density, Dirichlet boundaries") {
    CylindricalPoissonSolver2D::DomainProp prop;
    prop.extents = IntVec<2>{3, 3};
    prop.dx = Vec<2>{1.0, 1.0};

    std::vector<CylindricalPoissonSolver2D::Region> regions = {
        {CellType::BoundaryDirichlet, {0, 0}, {0, 2}, []{ return 1.0; }},
        {CellType::BoundaryDirichlet, {2, 0}, {2, 2}, []{ return 2.0; }},
        {CellType::BoundaryDirichlet, {0, 0}, {2, 0}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {0, 2}, {2, 2}, []{ return 3.0; }},
    };

    CylindricalPoissonSolver2D solver(prop, regions);

    Matrix<2> rho({3, 3});
    rho.fill(0.0);
    Matrix<2> phi;

    solver.solve(phi, rho);

    REQUIRE(phi(0, 1) == Catch::Approx(1.0)); 
    REQUIRE(phi(2, 1) == Catch::Approx(2.0));
    REQUIRE(phi(1, 0) == Catch::Approx(0.0));
    REQUIRE(phi(1, 2) == Catch::Approx(3.0));

    REQUIRE(phi(0, 0) == Catch::Approx(0.0));
    REQUIRE(phi(2, 0) == Catch::Approx(0.0));
    REQUIRE(phi(0, 2) == Catch::Approx(3.0));
    REQUIRE(phi(2, 2) == Catch::Approx(3.0));
}

TEST_CASE("CylindricalPoissonSolver2D: Constant density, Dirichlet boundaries") {
    CylindricalPoissonSolver2D::DomainProp prop;
    prop.extents = IntVec<2>{3, 3};
    prop.dx = Vec<2>{1.0, 1.0};

    std::vector<CylindricalPoissonSolver2D::Region> regions = {
        {CellType::BoundaryDirichlet, {0, 0}, {0, 2}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {2, 0}, {2, 2}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {0, 0}, {2, 0}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {0, 2}, {2, 2}, []{ return 0.0; }},
    };

    CylindricalPoissonSolver2D solver(prop, regions);

    Matrix<2> rho({3});
    rho.fill(1.0);
    Matrix<2> phi;

    solver.solve(phi, rho);

    double center = phi(1, 1);
    for (int i = 1; i < 2; ++i)
        for (int j = 1; j < 2; ++j)
            REQUIRE(phi(i, j) == Catch::Approx(center));
}

TEST_CASE("StructPoissonSolver2D: Analytical solution") {
    // D2(phi) = d2/dx2(phi) + d2/dy2(phi)
    // d2/dx2(phi) = -(pi/Lx)^2 * phi
    // d2/dy2(phi) = -(pi/Ly)^2 * phi
    // D2(phi) = -pi^2 * (1/Lx^2 + 1/Ly^2) * phi
    // rho = -eps0 * D2(phi) = eps0 * pi^2 * (1/Lx^2 + 1/Ly^2) * phi

    const int nx = 17, ny = 17;
    const double Lx = 1.0, Ly = 1.0;

    StructPoissonSolver2D::DomainProp prop;
    prop.extents = IntVec<2>{nx, ny};
    prop.dx = Vec<2>{Lx / (nx - 1), Ly / (ny - 1)};

    std::vector<StructPoissonSolver2D::Region> regions = {
        {CellType::BoundaryDirichlet, {0, 0}, {0, ny - 1}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {nx - 1, 0}, {nx - 1, ny - 1}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {0, 0}, {nx - 1, 0}, []{ return 0.0; }},
        {CellType::BoundaryDirichlet, {0, ny - 1}, {nx - 1, ny - 1}, []{ return 0.0; }},
    };

    StructPoissonSolver2D solver(prop, regions);

    Matrix<2> rho({(size_t)nx, (size_t)ny});
    const double A = spark::constants::eps0 * spark::constants::pi * spark::constants::pi * (1.0/(Lx*Lx) + 1.0/(Ly*Ly));
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = i * prop.dx.x;
            double y = j * prop.dx.y;
            rho(i, j) = A * sin(spark::constants::pi * x / Lx) * sin(spark::constants::pi * y / Ly);
        }
    }
    Matrix<2> phi_numerical;
    solver.solve(phi_numerical, rho);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = i * prop.dx.x;
            double y = j * prop.dx.y;
            double phi_analytical = sin(spark::constants::pi * x / Lx) * sin(spark::constants::pi * y / Ly);
            REQUIRE(phi_numerical(i, j) == Catch::Approx(phi_analytical).margin(1e-2));
        }
    }
}

TEST_CASE("CylindricalPoissonSolver2D: Analytical solution (r^2)") {
    // D2(phi) = d2/dr2(phi) + (1/r)*d/dr(phi) + d2/dz2(phi)
    // phi = r^2, D2(phi) = 4.
    // rho = -eps0 * D2(phi)
    // rho = -4 * eps0

    const int nz = 11, nr = 11;
    const double Lz = 1.0, R = 1.0;

    CylindricalPoissonSolver2D::DomainProp prop;
    prop.extents = IntVec<2>{nz, nr};
    prop.dx = Vec<2>{Lz / (nz - 1), R / (nr - 1)};

    std::vector<CylindricalPoissonSolver2D::Region> regions;

    for(int j=0; j < nr; ++j) {
        double r = j * prop.dx.y;
        regions.push_back({CellType::BoundaryDirichlet, {0, j}, {0, j}, [r]{ return r*r; }});
    }
    for(int j=0; j < nr; ++j) {
        double r = j * prop.dx.y;
        regions.push_back({CellType::BoundaryDirichlet, {nz-1, j}, {nz-1, j}, [r]{ return r*r; }});
    }
    regions.push_back({CellType::BoundaryDirichlet, {0, nr - 1}, {nz - 1, nr - 1}, [R]{ return R*R; }});

    CylindricalPoissonSolver2D solver(prop, regions);

    Matrix<2> rho({(size_t)nz, (size_t)nr});
    rho.fill(-4.0 * spark::constants::eps0);
    
    Matrix<2> phi_numerical;
    solver.solve(phi_numerical, rho);

    for (int i = 0; i < nz; ++i) {
        for (int j = 0; j < nr; ++j) {
            double r_coord = j * prop.dx.y;
            double phi_analytical = r_coord * r_coord;
            REQUIRE(phi_numerical(i, j) == Catch::Approx(phi_analytical).margin(1e-3));
        }
    }
}