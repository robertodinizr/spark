#pragma once

#include <cmath>

namespace spark::core {

template <typename T, unsigned N>
struct TVec {
    T norm() const;
    TVec normalized() const;
};

template <typename T>
struct TVec<T, 1> {
    T x{};
    T norm() const { return x; }
    T sum() const { return x; }
    T mul() const { return x; }
    TVec normalized() const { return TVec{x / norm()}; }

    template <typename G>
    TVec<G, 1> to() const {
        return TVec<G, 1>{static_cast<G>(x)};
    }
};

template <typename T>
struct TVec<T, 2> {
    T x{}, y{};
    T norm() const { return std::sqrt(x * x + y * y); }
    T sum() const { return x + y; }
    T mul() const { return x * y; }
    TVec normalized() const {
        const T n = norm();
        return TVec{x / n, y / n};
    }

    template <typename G>
    TVec<G, 2> to() const {
        return TVec<G, 2>{static_cast<G>(x), static_cast<G>(y)};
    }
};

template <typename T>
struct TVec<T, 3> {
    T x{}, y{}, z{};
    T norm() const { return std::sqrt(x * x + y * y + z * z); }
    T sum() const { return x + y + z; }
    T mul() const { return x * y * z; }
    TVec normalized() const {
        const T n = norm();
        return TVec{x / n, y / n, z / n};
    }

    template <typename G>
    TVec<G, 3> to() const {
        return TVec<G, 3>{static_cast<G>(x), static_cast<G>(y), static_cast<G>(z)};
    }
};

template <unsigned N>
using Vec = TVec<double, N>;

template <unsigned N>
using IntVec = TVec<int, N>;

template <unsigned N>
using ULongVec = TVec<size_t, N>;

}  // namespace spark::core