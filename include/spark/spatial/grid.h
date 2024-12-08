#pragma once

#include <spark/core/vec.h>

#include <vector>

namespace spark::spatial {

template <unsigned N>
class UniformGrid {
public:
    UniformGrid() = default;
    UniformGrid(const std::array<double, N>& l, const std::array<size_t, N>& n) {
        l_ = l;
        n_ = n;
        n_total_ = 1;
        for (unsigned i = 0; i < N; i++) {
            dx_[i] = l[i] / static_cast<double>(n[i] - 1);
            n_total_ *= n_[i];
        }

        data_.resize(n_total_);
        set(0.0);
    }

    void set(double v) { std::fill(data_.begin(), data_.end(), v); }
    std::vector<double>& data() { return data_; }
    const std::vector<double>& data() const { return data_; }
    double* data_ptr() const { return const_cast<double*>(data_.data()); }

    auto n_total() const { return n_total_; }
    auto n() const { return n_; }
    auto l() const { return l_; }
    auto dx() const { return dx_; }

    void apply(double mul, double add) {
        for (size_t i = 0; i < n_total_; i++) {
            data_[i] = data_[i] * mul + add;
        }
    }

private:
    size_t n_total_ = 0;
    core::ULongVec<N> n_;
    core::Vec<N> l_, dx_;
    std::vector<double> data_;
};

class AverageGrid {
public:
    AverageGrid() = default;
    AverageGrid(const UniformGrid<1>& grid) : m_average_grid(grid) { m_average_grid.set(0); }

    void add(const UniformGrid<1>& grid) {
        auto& av = m_average_grid.data();
        auto& gr = grid.data();

        for (size_t i = 0; i < m_average_grid.n_total(); i++) {
            av[i] =
                (av[i] * (double)m_count / (double)(m_count + 1)) + (gr[i] / (double)(m_count + 1));
        }

        m_count++;
    }
    std::vector<double>& get() { return m_average_grid.data(); }

private:
    UniformGrid<1> m_average_grid;
    size_t m_count = 0;
};

}  // namespace spark::spatial