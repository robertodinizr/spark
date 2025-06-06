#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <limits>

#include "spark/core/vec.h"

using namespace spark::core;

static double square(double x) {
    return x * x;
}

TEST_CASE("TVec<T, 1>") {
    SECTION("Construction and initialization") {
        Vec<1> v;
        REQUIRE(v.x == 0.0);

        Vec<1> v_init{10.5};
        REQUIRE(v_init.x == Catch::Approx(10.5));
    }

    SECTION("Vector properties") {
        Vec<1> v{-5.0};
        REQUIRE(v.norm() == Catch::Approx(5.0));
        REQUIRE(v.sum() == Catch::Approx(-5.0));
        REQUIRE(v.mul() == Catch::Approx(-5.0));

        Vec<1> v_norm = v.normalized();
        REQUIRE(v_norm.x == Catch::Approx(-1.0));
    }

    SECTION("Arithmetic operators") {
        Vec<1> a{10.0};
        Vec<1> b{2.0};

        REQUIRE((a + b).x == Catch::Approx(12.0));
        REQUIRE((a - b).x == Catch::Approx(8.0));
        REQUIRE((a * b).x == Catch::Approx(20.0));
        REQUIRE((a / b).x == Catch::Approx(5.0));
    }
}


TEST_CASE("TVec<T, 2>") {
    
    SECTION("Construction and initialization") {
        Vec<2> v;
        REQUIRE(v.x == 0.0);
        REQUIRE(v.y == 0.0);

        Vec<2> v_init{3.0, -4.0};
        REQUIRE(v_init.x == Catch::Approx(3.0));
        REQUIRE(v_init.y == Catch::Approx(-4.0));
    }

    SECTION("Vector properties") {
        Vec<2> v{3.0, 4.0};
        REQUIRE(v.norm() == Catch::Approx(5.0));
        REQUIRE(v.sum() == Catch::Approx(7.0));
        REQUIRE(v.mul() == Catch::Approx(12.0));
    }

    SECTION("Normalization") {
        Vec<2> v{3.0, -4.0};
        Vec<2> v_norm = v.normalized();
        
        REQUIRE(v_norm.x == Catch::Approx(0.6));
        REQUIRE(v_norm.y == Catch::Approx(-0.8));
        REQUIRE(v_norm.norm() == Catch::Approx(1.0));

        Vec<2> zero_vec{0.0, 0.0};
        Vec<2> zero_norm = zero_vec.normalized();
        REQUIRE(std::isnan(zero_norm.x));
        REQUIRE(std::isnan(zero_norm.y));
    }

    SECTION("Arithmetic Operators") {
        Vec<2> a{10.0, 20.0};
        Vec<2> b{2.0, 5.0};

        Vec<2> sum = a + b;
        REQUIRE(sum.x == Catch::Approx(12.0));
        REQUIRE(sum.y == Catch::Approx(25.0));
        
        Vec<2> sub = a - b;
        REQUIRE(sub.x == Catch::Approx(8.0));
        REQUIRE(sub.y == Catch::Approx(15.0));
        
        Vec<2> mul = a * b;
        REQUIRE(mul.x == Catch::Approx(20.0));
        REQUIRE(mul.y == Catch::Approx(100.0));
        
        Vec<2> div = a / b;
        REQUIRE(div.x == Catch::Approx(5.0));
        REQUIRE(div.y == Catch::Approx(4.0));
        
        Vec<2> scaled = a * 2.0;
        REQUIRE(scaled.x == Catch::Approx(20.0));
        REQUIRE(scaled.y == Catch::Approx(40.0));
        
        Vec<2> scaled_rev = 2.0 * a;
        REQUIRE(scaled.x == Catch::Approx(scaled_rev.x));
        REQUIRE(scaled.y == Catch::Approx(scaled_rev.y));
    }

    SECTION("Comparison") {
        Vec<2> a{1.0, 2.0};
        Vec<2> b{1.0, 2.0};
        Vec<2> c{1.0, 3.0};

        REQUIRE(a == b);
        REQUIRE(a != c);
    }

    SECTION("Type conversion") {
        Vec<2> v_double{10.7, -3.9};
        IntVec<2> v_int = v_double.to<int>();

        REQUIRE(v_int.x == 10);
        REQUIRE(v_int.y == -3);
    }

    SECTION("Array conversion") {
        Vec<2> v{5.5, 6.6};
        auto v_arr = v.arr();
        REQUIRE(v_arr[0] == Catch::Approx(5.5));
        REQUIRE(v_arr[1] == Catch::Approx(6.6));
    }

    SECTION("Apply method") {
        Vec<2> v{2.0, 3.0};
        Vec<2> v_sq = v.apply<square>();
        REQUIRE(v_sq.x == Catch::Approx(4.0));
        REQUIRE(v_sq.y == Catch::Approx(9.0));
    }
}


TEST_CASE("TVec<T, 3>") {
    
    SECTION("Construction and initialization") {
        Vec<3> v;
        REQUIRE(v.x == 0.0);
        REQUIRE(v.y == 0.0);
        REQUIRE(v.z == 0.0);

        Vec<3> v_init{1.0, 2.0, -3.0};
        REQUIRE(v_init.x == Catch::Approx(1.0));
        REQUIRE(v_init.y == Catch::Approx(2.0));
        REQUIRE(v_init.z == Catch::Approx(-3.0));
    }

    SECTION("Vector properties") {
        Vec<3> v{2.0, 3.0, 4.0};
        REQUIRE(v.norm() == Catch::Approx(std::sqrt(29.0)));
        REQUIRE(v.sum() == Catch::Approx(9.0));
        REQUIRE(v.mul() == Catch::Approx(24.0));
    }

    SECTION("Normalization") {
        Vec<3> v{2.0, -3.0, 6.0};
        REQUIRE(v.norm() == Catch::Approx(7.0));

        Vec<3> v_norm = v.normalized();
        REQUIRE(v_norm.x == Catch::Approx(2.0 / 7.0));
        REQUIRE(v_norm.y == Catch::Approx(-3.0 / 7.0));
        REQUIRE(v_norm.z == Catch::Approx(6.0 / 7.0));
        REQUIRE(v_norm.norm() == Catch::Approx(1.0));
    }

    SECTION("Arithmetic operators") {
        Vec<3> a{10.0, 20.0, 30.0};
        Vec<3> b{2.0, 4.0, 5.0};

        Vec<3> sum = a + b;
        REQUIRE(sum.x == Catch::Approx(12.0));
        REQUIRE(sum.y == Catch::Approx(24.0));
        REQUIRE(sum.z == Catch::Approx(35.0));

        Vec<3> scaled = a / 2.0;
        REQUIRE(scaled.x == Catch::Approx(5.0));
        REQUIRE(scaled.y == Catch::Approx(10.0));
        REQUIRE(scaled.z == Catch::Approx(15.0));
    }
}