#pragma once

#include <vector>
#include <functional>

#include "kn/core/vec.h"

namespace kn::particle {

    class ChargedSpecies1D1V {
    public:
        ChargedSpecies1D1V() = default;
        ChargedSpecies1D1V(double q, double m) : m_q(q), m_m(m) {}
        
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
        
        ChargedSpecies1D3V() = default;
        ChargedSpecies1D3V(double q, double m) : m_q(q), m_m(m) {}
        
        void add(size_t n);
        void add(size_t n, std::function<void(core::Vec3&, double&)> sampler);
        void add_copy(size_t idx);
        void remove(size_t idx);

        double m() const { return m_m; }
        core::Vec3* v() const { return (core::Vec3*) m_v.data(); };
        double* x() const { return (double*) m_x.data(); };
        double* f() const { return (double*) m_f.data(); };
        size_t n() const { return m_n; }
        double q() const { return m_q; }

    private:
        std::vector<core::Vec3> m_v;
        std::vector<double> m_x;
        std::vector<double> m_f;
        size_t m_n = 0;

        double m_q = 0;
        double m_m = 0;
    };

    template <unsigned NX, unsigned NV>
    class ChargedSpecies {
    public:
        ChargedSpecies() = default;
        ChargedSpecies(double q, double m) : m_q(q), m_m(m) {}
        
        void add(size_t n);
        void add(size_t n, std::function<void(core::Vec<NV>&, core::Vec<NX>&)> sampler);
        void add_copy(size_t idx);
        void remove(size_t idx);
    
        double m() const { return m_m; }
        core::Vec<NV>* v() const { return (core::Vec<NV>*) m_v.data(); };
        core::Vec<NX>* x() const { return (core::Vec<NX>*) m_x.data(); };
        core::Vec<NX>* f() const { return (core::Vec<NX>*) m_f.data(); };
        size_t n() const { return m_n; }
        double q() const { return m_q; }

    private:
        std::vector<core::Vec<NV>> m_v;
        std::vector<core::Vec<NX>> m_x;
        std::vector<core::Vec<NX>> m_f;
        size_t m_n = 0;

        double m_q = 0;
        double m_m = 0;
    };
}