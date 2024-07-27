#pragma once

#include <vector>

namespace kn::spatial {

    class UniformGrid {
    public:
        UniformGrid() = default;
        UniformGrid(double l, size_t n);

        void set(double v);
        std::vector<double>& data();
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
}