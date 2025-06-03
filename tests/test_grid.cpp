#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include "spark/spatial/grid.h"
#include "spark/core/vec.h"

using namespace spark;
using namespace spark::spatial;

TEST_CASE("TUniformGrid basic construction and properties") {
    core::Vec<2> l{10.0, 5.0};
    core::ULongVec<2> n{11, 6};
    TUniformGrid<double, 2> grid(l, n);

    REQUIRE(grid.n().x == 11);
    REQUIRE(grid.n().y == 6);
    REQUIRE(grid.l().x == Catch::Approx(10.0));
    REQUIRE(grid.l().y == Catch::Approx(5.0));
    REQUIRE(grid.dx().x == Catch::Approx(10.0 / 10.0));
    REQUIRE(grid.dx().y == Catch::Approx(5.0 / 5.0));
    REQUIRE(grid.n_total() == 66);
}

TEST_CASE("TUniformGrid set and data access") {
    core::Vec<2> l{2.0, 2.0};
    core::ULongVec<2> n{3, 3};
    TUniformGrid<double, 2> grid(l, n);

    grid.set(42.0);
    auto& data = grid.data();
    for (size_t i = 0; i < data.count(); ++i) {
        REQUIRE(data[i] == Catch::Approx(42.0));
    }
}

TEST_CASE("TUniformGrid apply mul/add") {
    core::Vec<2> l{2.0, 2.0};
    core::ULongVec<2> n{2, 2};
    TUniformGrid<double, 2> grid(l, n);

    grid.set(2.0);
    grid.apply(3.0, 1.0); // (2*3)+1 = 7
    auto& data = grid.data();
    for (size_t i = 0; i < data.count(); ++i) {
        REQUIRE(data[i] == Catch::Approx(7.0));
    }
}

TEST_CASE("TUniformGrid radius and volume calculations") {
    core::Vec<2> l{4.0, 4.0};
    core::ULongVec<2> n{3, 3};
    TUniformGrid<double, 2> grid(l, n);

    double r0 = grid.get_radius(0);
    double r1 = grid.get_radius(1);
    double mid_r0 = grid.get_mid_radius(0);

    REQUIRE(r0 == Catch::Approx(0.0));
    REQUIRE(r1 == Catch::Approx(grid.dx().y));
    REQUIRE(mid_r0 == Catch::Approx(0.5 * grid.dx().y));

    double vol = grid.get_annular_cell_volume(1, 1);
    REQUIRE(vol > 0.0);
}

TEST_CASE("TUniformGrid get_cell_area") {
    core::Vec<2> l{4.0, 4.0};
    core::ULongVec<2> n{3, 3};
    TUniformGrid<double, 2> grid(l, n);

    double area = 0.0;
    grid.get_cell_area(1, 1, area);
    REQUIRE(area > 0.0);
}

TEST_CASE("AverageGrid add and get") {
    core::Vec<2> l{2.0, 2.0};
    core::ULongVec<2> n{2, 2};
    UniformGrid<2> grid1(l, n);
    UniformGrid<2> grid2(l, n);

    grid1.set(2.0);
    grid2.set(4.0);

    AverageGrid<2> avg(grid1);
    avg.add(grid1); // count = 1, avg = 2
    avg.add(grid2); // count = 2, avg = (2+4)/2 = 3

    auto& data = avg.get();
    for (size_t i = 0; i < data.count(); ++i) {
        REQUIRE(data[i] == Catch::Approx(3.0));
    }
}