#pragma once

#include <vector>

namespace spark::spatial {

class UniformGrid {
public:
    UniformGrid() = default;
    UniformGrid(double l, size_t n);

    void set(double v);
    std::vector<double>& data();
    const std::vector<double>& data() const;
    double* data_ptr() const;

    size_t n() const;
    double l() const;
    double dx() const;

    void apply(double mul, double add);

private:
    size_t m_n = 0;
    double m_l = 0.0, m_dx = 0.0;
    std::vector<double> m_data;
};

class AverageGrid {
public:
    AverageGrid() = default;
    AverageGrid(double l, size_t n) : m_average_grid(l, n) {}
    AverageGrid(const UniformGrid& grid) : m_average_grid(grid.l(), grid.n()) {}

    void add(const UniformGrid& grid);
    std::vector<double>& get() { return m_average_grid.data(); }

private:
    UniformGrid m_average_grid;
    size_t m_count = 0;
};

}  // namespace spark::spatial