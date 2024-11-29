#pragma once

#include <functional>
#include <vector>

#include "spark/core/vec.h"

namespace spark::particle {

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
    core::Vec<NV>* v() const { return (core::Vec<NV>*)m_v.data(); };
    core::Vec<NX>* x() const { return (core::Vec<NX>*)m_x.data(); };
    core::Vec<NX>* f() const { return (core::Vec<NX>*)m_f.data(); };
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

}  // namespace spark::particle