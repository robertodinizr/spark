#pragma once

#include <vector>

#include "spark/core/matrix.h"
#include "spark/core/vec.h"

namespace spark::spatial {

template <unsigned N>
struct GridProp {
    core::Vec<N> l_;
    core::Vec<N> dx_;
};

template <unsigned N>
class UniformGrid {
public:
    UniformGrid() = default;
    UniformGrid(const core::Vec<N>& l, const core::ULongVec<N>& n) {
        l_ = l;
        data_.resize(n);
        dx_ = l / (n.template to<double>() - 1.0);
        set(0.0);
    }

    void set(double v) { data_.fill(v); }
    auto& data() { return data_; }
    const auto& data() const { return data_; }
    double* data_ptr() const { return const_cast<double*>(data_.data_ptr()); }

    auto n_total() const { return data_.count(); }
    auto n() const { return data_.size(); }
    auto l() const { return l_; }
    auto dx() const { return dx_; }

    void apply(double mul, double add) {
        const auto count = data_.count();
        for (size_t i = 0; i < count; i++) {
            data_[i] = data_[i] * mul + add;
        }
    }

private:
    core::Vec<N> l_, dx_;
    core::Matrix<N> data_;
};

template <unsigned N>
class AverageGrid {
public:
    AverageGrid() = default;
    AverageGrid(const UniformGrid<N>& grid) : m_average_grid(grid) { m_average_grid.set(0); }

    void add(const UniformGrid<1>& grid) {
        auto& av = m_average_grid.data();
        auto& gr = grid.data();

        const auto count = m_average_grid.n_total();
        for (size_t i = 0; i < count; i++) {
            av(i) = (av(i) * static_cast<double>(m_count) / static_cast<double>(m_count + 1)) +
                    (gr(i) / static_cast<double>(m_count + 1));
        }

        m_count++;
    }

    auto& get() { return m_average_grid.data(); }

private:
    UniformGrid<N> m_average_grid;
    size_t m_count = 0;
};

}  // namespace spark::spatial