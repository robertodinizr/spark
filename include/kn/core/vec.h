#pragma once

#include <cmath>

namespace kn::core {

    struct Vec2 {
        double x = 0.0, y = 0.0;
        
        double norm() const {
            return std::sqrt(x * x + y * y);
        }

        Vec2 normalized() const {
            double n = norm();
            return Vec2{x / n, y / n};
        }
    };

    struct Vec3 {
        double x = 0.0, y = 0.0, z = 0.0;
        
        double norm() const {
            return std::sqrt(x * x + y * y + z * z);
        }

        Vec3 normalized() const {
            double n = norm();
            return Vec3{x / n, y / n, z / n};
        }
    };
}