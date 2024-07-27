#pragma once

#include <vector>
#include <functional>

namespace kn::particle {

    class ChargedSpecies {
    public:
        ChargedSpecies() = default;
        ChargedSpecies(double q, double m) : m_q(q), m_m(m) {}
        
        void add(size_t n);
        void add(size_t n, std::function<void(double&, double&)> sampler);
        
        double* v() const { return (double*) m_v.data(); };
        double* x() const { return (double*) m_x.data(); };
        double* f() const { return (double*) m_f.data(); };
        size_t n() const { return m_n; }

    private:
        std::vector<double> m_v;
        std::vector<double> m_x;
        std::vector<double> m_f;
        size_t m_n = 0;

        double m_q = 0;
        double m_m = 0;
    };
}