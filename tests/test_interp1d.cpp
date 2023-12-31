#include "sci.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <vector>

TEST_CASE( "Nearest interpolation", "[interp1d]" ) {
    
    auto x = sci::buffer<double>({1.0, 2.0, 3.0, 4.0});
    auto y = sci::buffer<double>({1.0, 2.0, 3.0, 4.0});

    auto nearest = sci::interpolate::interp1d<double>(
        std::move(x),
        std::move(y),
        sci::interpolate::interp1d_type::Nearest
    );

    REQUIRE(nearest(2.5) == 2.0);
    REQUIRE(isnan(nearest(-10.0)));
    REQUIRE(nearest(5.0) == 4.0);
    REQUIRE(nearest(3.999) == 3.0);
    REQUIRE(nearest(1.0) == 1.0);
    REQUIRE(nearest(4.0) == 4.0);


    sci::buffer<int> x2 = {2, 4, 8, 16, 32};
    sci::buffer<int> y2 = {2, 4, 8, 16, 32};

    auto nearest_int = sci::interpolate::interp1d<int>(
        std::move(x2), std::move(y2), 
        sci::interpolate::interp1d_type::Nearest
    );

    REQUIRE(nearest_int(3) == 2);

}

TEST_CASE("Nearest up interpoation", "[interp1d]")
{
    auto x = sci::buffer<double>({1.0, 2.0, 3.0, 4.0});
    auto y = sci::buffer<double>({1.0, 2.0, 3.0, 4.0});

    auto nearest = sci::interpolate::interp1d<double>(
        std::move(x),
        std::move(y),
        sci::interpolate::interp1d_type::NearestUp
    );

    REQUIRE(nearest(2.5) == 3.0);
    REQUIRE(isnan(nearest(5.0)));
    REQUIRE(nearest(3.999) == 4.0);
    REQUIRE(nearest(1.0) == 1.0);
    REQUIRE(nearest(4.0) == 4.0); 
}


