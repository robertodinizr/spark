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
        
        double m() const { return m_m; }
        double* v() const { return (double*) m_v.data(); };
        double* x() const { return (double*) m_x.data(); };
        double* f() const { return (double*) m_f.data(); };
        size_t n() const { return m_n; }
        size_t q() const { return m_q; }

    private:
        std::vector<double> m_v;
        std::vector<double> m_x;
        std::vector<double> m_f;
        size_t m_n = 0;

        double m_q = 0;
        double m_m = 0;
    };

    class ChargedSpecies1D3V {
    public:
        struct Vec3 {
            double x = 0.0, y = 0.0, z = 0.0;
            
            double norm() {
                return std::sqrt(x * x + y * y + z * z);
            }

            Vec3 normalized() {
                double n = norm();
                return Vec3{x / n, y / n, z / n};
            }
        };

        ChargedSpecies1D3V() = default;
        ChargedSpecies1D3V(double q, double m) : m_q(q), m_m(m) {}
        
        void add(size_t n);
        void add(size_t n, std::function<void(Vec3&, double&)> sampler);
        void remove(size_t idx);

        double m() const { return m_m; }
        Vec3* v() const { return (Vec3*) m_v.data(); };
        double* x() const { return (double*) m_x.data(); };
        double* f() const { return (double*) m_f.data(); };
        size_t n() const { return m_n; }
        size_t q() const { return m_q; }

    private:
        std::vector<Vec3> m_v;
        std::vector<double> m_x;
        std::vector<double> m_f;
        size_t m_n = 0;

        double m_q = 0;
        double m_m = 0;
    };
}