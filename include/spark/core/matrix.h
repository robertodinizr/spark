
#ifndef SPARK_MATRIX_H
#define SPARK_MATRIX_H

#include <vector>

#include "vec.h"

namespace spark::core {

#define ENABLE_IF_(NN) template <unsigned N_ = N, std::enable_if_t<N_ == NN>* = nullptr>

template <typename T, unsigned N>
    requires(N >= 1 && N <= 3)
class TMatrix {
public:
    TMatrix() = default;
    explicit TMatrix(const ULongVec<N>& size) { resize(size); }

    void resize(const ULongVec<N>& size) {
        size_ = size;
        count_ = size_.mul();
        data_.resize(count_);
    }

    void fill(const T& value) { std::fill(data_.begin(), data_.end(), value); }

    ENABLE_IF_(2)
    void fill(const T& value, const ULongVec<N>& lower_left, const ULongVec<N>& upper_right) {
        for (size_t i = 0; i < upper_right.x - lower_left.x + 1; ++i)
            for (size_t j = 0; j < upper_right.y - lower_left.y + 1; ++j)
                operator()(i + lower_left.x, j + lower_left.y) = value;
    }

    ENABLE_IF_(1)
    void fill(const T& value, const ULongVec<N>& x0, const ULongVec<N>& x1) {
        for (size_t i = x0.x; i < x1.x - x0.x + 1; ++i)
            data_[i] = value;
    }

    std::vector<T>& data() { return data_; }
    const std::vector<T>& data() const { return data_; }
    T* data_ptr() const { return const_cast<double*>(data_.data()); }

    auto begin() { return data_.begin(); }
    auto end() { return data_.end(); }

    ULongVec<N> size() const { return size_; }
    size_t count() const { return count_; }

#define IDX2D(i, j) ((i) * size_.y + (j))
#define IDX3D(i, j, k) ((i) * size_.y * size_.z + (j) * size_.z + (k))

    ENABLE_IF_(1) size_t index(const size_t i) const { return i; }
    ENABLE_IF_(2) size_t index(const size_t i, const size_t j) const { return IDX2D(i, j); }
    ENABLE_IF_(3) size_t index(const size_t i, const size_t j, const size_t k) const {
        return IDX3D(i, j, k);
    }

    ENABLE_IF_(1) T operator()(const size_t i) const { return data_[i]; }
    ENABLE_IF_(1) T& operator()(const size_t i) { return data_[i]; }
    ENABLE_IF_(2) T operator()(const size_t i, const size_t j) const { return data_[IDX2D(i, j)]; }
    ENABLE_IF_(2) T& operator()(const size_t i, const size_t j) { return data_[IDX2D(i, j)]; }
    ENABLE_IF_(3) T operator()(const size_t i, const size_t j, const size_t k) const {
        return data_[IDX3D(i, j, k)];
    }
    ENABLE_IF_(3) T& operator()(const size_t i, const size_t j, const size_t k) {
        return data_[IDX3D(i, j, k)];
    }

    T operator[](const ULongVec<N>& ij) const {
        if constexpr (N == 1)
            return ij.x;
        else if constexpr (N == 2)
            return IDX2D(ij.x, ij.y);
        else
            return IDX3D(ij.x, ij.y, ij.k);
    }

    T& operator[](const ULongVec<N>& ij) {
        if constexpr (N == 1)
            return ij.x;
        else if constexpr (N == 2)
            return IDX2D(ij.x, ij.y);
        else
            return IDX3D(ij.x, ij.y, ij.k);
    }

#undef IDX2D
#undef IDX3D

private:
    size_t count_ = 0;
    ULongVec<N> size_;
    std::vector<T> data_;
};

#undef ENABLE_IF_

template <unsigned N>
using Matrix = TMatrix<double, N>;

}  // namespace spark::core

#endif  // MATRIX_H
