// tests/test_matrix.cpp
#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "spark/core/matrix.h"
#include "spark/core/vec.h"

using namespace spark::core;

TEST_CASE("TMatrix<T, 1>") {
    SECTION("Construction and sizing") {
        Matrix<1> m_empty;
        REQUIRE(m_empty.size().x == 0);
        REQUIRE(m_empty.count() == 0);

        Matrix<1> m({10});
        REQUIRE(m.size().x == 10);
        REQUIRE(m.count() == 10);

        m.resize({20});
        REQUIRE(m.size().x == 20);
        REQUIRE(m.count() == 20);
    }

    SECTION("Data access") {
        Matrix<1> m({5});
        m.fill(3.14);

        for (size_t i = 0; i < m.count(); ++i) {
            REQUIRE(m(i) == Catch::Approx(3.14));
            REQUIRE(m[i] == Catch::Approx(3.14));
        }

        m(2) = 1.61;
        REQUIRE(m(2) == Catch::Approx(1.61));

        auto* ptr = m.data_ptr();
        REQUIRE(ptr[2] == Catch::Approx(1.61));
    }
}

TEST_CASE("TMatrix<T, 2>") {
    SECTION("Construction and sizing") {
        Matrix<2> m({10, 5});
        REQUIRE(m.size().x == 10);
        REQUIRE(m.size().y == 5);
        REQUIRE(m.count() == 50);

        m.resize({20, 10});
        REQUIRE(m.size().x == 20);
        REQUIRE(m.size().y == 10);
        REQUIRE(m.count() == 200);
    }

    SECTION("Data access") {
        Matrix<2> m({4, 3});
        m.fill(1.0);

        m(2, 1) = 5.5;
        REQUIRE(m(2, 1) == Catch::Approx(5.5));
        ULongVec<2> coords = {3, 2};
        m[coords] = 6.6;
        REQUIRE(m[coords] == Catch::Approx(6.6));
        REQUIRE(m(3, 2) == Catch::Approx(6.6));
    }

    SECTION("Matrix fill") {
        Matrix<2> m({10, 10});
        m.fill(0.0);

        IntVec<2> lower_left = {2, 2};
        IntVec<2> upper_right = {5, 5};
        m.fill(9.9, lower_left, upper_right);

        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                if (i >= 2 && i <= 5 && j >= 2 && j <= 5) {
                    REQUIRE(m(i, j) == Catch::Approx(9.9));
                } else {
                    REQUIRE(m(i, j) == Catch::Approx(0.0));
                }
            }
        }
    }
}

TEST_CASE("TMatrix<T, 3>") {
    SECTION("Construction and sizing") {
        Matrix<3> m({10, 5, 2});
        REQUIRE(m.size().x == 10);
        REQUIRE(m.size().y == 5);
        REQUIRE(m.size().z == 2);
        REQUIRE(m.count() == 100);

        m.resize({3, 3, 3});
        REQUIRE(m.size().x == 3);
        REQUIRE(m.size().y == 3);
        REQUIRE(m.size().z == 3);
        REQUIRE(m.count() == 27);
    }

    SECTION("Data access") {
        Matrix<3> m({2, 2, 2});
        m.fill(1.0);

        m(1, 0, 1) = 8.8;
        REQUIRE(m(1, 0, 1) == Catch::Approx(8.8));
        REQUIRE(m(0, 0, 0) == Catch::Approx(1.0));
    }
}