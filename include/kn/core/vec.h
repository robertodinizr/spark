#pragma once

#include <cmath>

namespace kn::core {

    template <unsigned N> 
    struct Vec {
        double norm() const;
        Vec normalized() const;
    };

    template <>
    struct Vec<1> {
        double x = 0.0;
        double norm() const { return x; }
        Vec normalized() const { return Vec{x / norm()}; }
    };

    template <>
    struct Vec<2> {
        double x = 0.0, y = 0.0;
        double norm() const { return std::sqrt(x*x + y*y); }
        Vec normalized() const { 
            const double n = norm();
            return Vec{x / n, y / n}; 
        }
    };

    template <>
    struct Vec<3> {
        double x = 0.0, y = 0.0, z = 0.0;
        double norm() const { return std::sqrt(x*x + y*y + z*z); }
        Vec normalized() const { 
            const double n = norm();
            return Vec{x / n, y / n, z / n}; 
        }
    };

}