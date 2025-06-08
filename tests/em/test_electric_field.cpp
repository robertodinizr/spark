#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "spark/em/electric_field.h"
#include "spark/spatial/grid.h"
#include "spark/core/vec.h"

using namespace spark;
using namespace spark::em;
using namespace spark::spatial;
using namespace spark::core;

TEST_CASE("Electric field calculation (1D)") {
    // phi(x) = A*x + B, E = -d(phi)/dx = -A
    // A = 10.0, E = -10.0
    const double A = 10.0;
    const ULongVec<1> n = {11};
    const Vec<1> l = {1.0};
    UniformGrid<1> phi(l, n);
    
    for (size_t i = 0; i < n.x; ++i) {
        phi.data()[i] = A * static_cast<double>(i) * phi.dx().x;
    }

    TMatrix<Vec<1>, 1> e_field;
    electric_field<1>(phi, e_field);

    REQUIRE(e_field.count() == n.x);
    for (size_t i = 0; i < e_field.count(); ++i) {
        REQUIRE(e_field[i].x == Catch::Approx(-A));
    }
}

TEST_CASE("Electric field calculation (2D Cartesian)") {
    // phi(x, y) = Ax + By + C, Ex = -A, Ey = -B
    // A = 5.0, B = -3.0
    const double A = 5.0;
    const double B = -3.0;
    const ULongVec<2> n = {11, 11};
    const Vec<2> l = {1.0, 1.0};
    UniformGrid<2> phi(l, n);
    const auto dx = phi.dx();

    for (size_t i = 0; i < n.x; ++i) {
        for (size_t j = 0; j < n.y; ++j) {
            phi.data()(i, j) = A * static_cast<double>(i) * dx.x + B * static_cast<double>(j) * dx.y;
        }
    }

    TMatrix<Vec<2>, 2> e_field;
    electric_field<2>(phi, e_field);

    REQUIRE(e_field.size().x == n.x);
    REQUIRE(e_field.size().y == n.y);
    for (size_t i = 0; i < n.x; ++i) {
        for (size_t j = 0; j < n.y; ++j) {
            REQUIRE(e_field(i, j).x == Catch::Approx(-A));
            REQUIRE(e_field(i, j).y == Catch::Approx(-B));
        }
    }
}

TEST_CASE("Electric field calculation (2D Cylindrical)") {
    const ULongVec<2> n = {11, 12};
    const Vec<2> l = {1.0, 1.0};
    UniformGrid<2> phi(l, n);
    const auto dz = phi.dx().x;
    const auto dr = phi.dx().y;
    auto& phi_data = phi.data();

    TMatrix<Vec<2>, 2> e_field;

    SECTION("phi(z, r) = A*r^2") {
        // phi(z, r) = A*r^2 + C
        // Er = -d(phi)/dr = -2*A*r
        // Ez = -d(phi)/dz = 0
        const double A = 15.0;
        for (size_t i = 0; i < n.x; ++i) {
            for (size_t j = 0; j < n.y; ++j) {
                const double r = static_cast<double>(j) * dr;
                phi_data(i, j) = A * r * r;
            }
        }

        electric_field_cylindrical(phi, e_field);

        for (size_t i = 0; i < n.x; ++i) {
            for (size_t j = 0; j < n.y; ++j) {
                const double r = static_cast<double>(j) * dr;
                const double expected_Er = -2.0 * A * r;
                if (j == 0) {
                     REQUIRE(e_field(i, j).x == Catch::Approx(0.0));
                } else {
                     REQUIRE(e_field(i, j).x == Catch::Approx(expected_Er).margin(1e-9));
                }
                REQUIRE(e_field(i, j).y == Catch::Approx(0.0).margin(1e-9));
            }
        }
    }

    SECTION("phi(z, r) = B*z^2") {
        // phi(z, r) = B*z^2 + C
        // Er = -d(phi)/dr = 0
        // Ez = -d(phi)/dz = -2*B*z
        const double B = 25.0;
        for (size_t i = 0; i < n.x; ++i) {
            for (size_t j = 0; j < n.y; ++j) {
                const double z = static_cast<double>(i) * dz;
                phi_data(i, j) = B * z * z;
            }
        }

        electric_field_cylindrical(phi, e_field);

        for (size_t i = 0; i < n.x; ++i) {
            for (size_t j = 0; j < n.y; ++j) {
                const double z = static_cast<double>(i) * dz;
                const double expected_Ez = -2.0 * B * z;
                REQUIRE(e_field(i, j).x == Catch::Approx(0.0).margin(1e-9));
                REQUIRE(e_field(i, j).y == Catch::Approx(expected_Ez).margin(1e-9));
            }
        }
    }
    
    SECTION("phi(z, r) = A*r^2 + B*z^2") {
        // Er = -2*A*r
        // Ez = -2*B*z
        const double A = 15.0;
        const double B = 25.0;
        for (size_t i = 0; i < n.x; ++i) {
            for (size_t j = 0; j < n.y; ++j) {
                const double r = static_cast<double>(j) * dr;
                const double z = static_cast<double>(i) * dz;
                phi_data(i, j) = A * r * r + B * z * z;
            }
        }

        electric_field_cylindrical(phi, e_field);
        
        for (size_t i = 0; i < n.x; ++i) {
            for (size_t j = 0; j < n.y; ++j) {
                const double r = static_cast<double>(j) * dr;
                const double z = static_cast<double>(i) * dz;
                const double expected_Er = -2.0 * A * r;
                const double expected_Ez = -2.0 * B * z;

                if (j == 0) {
                     REQUIRE(e_field(i, j).x == Catch::Approx(0.0));
                } else {
                     REQUIRE(e_field(i, j).x == Catch::Approx(expected_Er).margin(1e-9));
                }
                REQUIRE(e_field(i, j).y == Catch::Approx(expected_Ez).margin(1e-9));
            }
        }
    }
}

TEST_CASE("Cylindrical field calculation with small grids") {
    TMatrix<Vec<2>, 2> e_field;

    SECTION("nz = 2") {
        const ULongVec<2> n = {2, 5};
        UniformGrid<2> phi({1.0, 1.0}, n);
        const double dz = phi.dx().x;
        const double B = 10.0;
        for(size_t i=0; i<n.x; ++i) for(size_t j=0; j<n.y; ++j) phi.data()(i,j) = B * (static_cast<double>(i)*dz);
        
        electric_field_cylindrical(phi, e_field);
        
        const double expected_Ez = -B;
        REQUIRE(e_field(0,0).y == Catch::Approx(expected_Ez));
        REQUIRE(e_field(1,0).y == Catch::Approx(expected_Ez));
    }

    SECTION("nr = 2") {
        const ULongVec<2> n = {5, 2};
        UniformGrid<2> phi({1.0, 1.0}, n);
        const double dr = phi.dx().y;
        const double A = 10.0;
        for(size_t i=0; i<n.x; ++i) for(size_t j=0; j<n.y; ++j) phi.data()(i,j) = A * (static_cast<double>(j)*dr);
        
        electric_field_cylindrical(phi, e_field);
        
        REQUIRE(e_field(0,0).x == Catch::Approx(0.0));
        const double expected_Er = -A;
        REQUIRE(e_field(0,1).x == Catch::Approx(expected_Er));
    }
    
    SECTION("nz = 1 or nr = 1") {
        const ULongVec<2> n_z1 = {1, 5};
        UniformGrid<2> phi_z1({1.0, 1.0}, n_z1);
        electric_field_cylindrical(phi_z1, e_field);
        REQUIRE(e_field(0,0).y == Catch::Approx(0.0));

        const ULongVec<2> n_r1 = {5, 1};
        UniformGrid<2> phi_r1({1.0, 1.0}, n_r1);
        electric_field_cylindrical(phi_r1, e_field);
        REQUIRE(e_field(0,0).x == Catch::Approx(0.0));
    }
}