#include "kn/spatial/grid.h"
#include <algorithm>

using namespace kn::spatial;

UniformGrid::UniformGrid(double l, size_t n) : 
    m_l(l), m_n(n) {
            m_dx = l / static_cast<double>(n - 1);
            m_data.resize(n);
            set(0.0);
};

void UniformGrid::set(double v) {
    std::fill(m_data.begin(), m_data.end(), v);
}

std::vector<double>& UniformGrid::data() {
    return m_data;
}

const std::vector<double>& UniformGrid::data() const {
    return m_data;
}

double* UniformGrid::data_ptr() const {
    return (double*) m_data.data();
}

size_t UniformGrid::n() const {
    return m_n;
}

double UniformGrid::l() const {
    return m_l;
}

double UniformGrid::dx() const {
    return m_dx;
}

void UniformGrid::apply(double mul, double add) {
    for(size_t i = 0; i < m_n; i++) {
        m_data[i] = m_data[i] * mul + add;
    }
}