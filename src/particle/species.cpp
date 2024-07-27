#include "kn/particle/species.h"

void kn::particle::ChargedSpecies::add(size_t n) {
    
    m_v.resize(m_n + n);
    m_x.resize(m_n + n);
    m_f.resize(m_n + n);
    
    for(size_t i = m_n; i < m_n + n; i++) {
        m_v[i] = 0;
        m_x[i] = 0;
        m_f[i] = 0;
    }

    m_n += n;
}

void kn::particle::ChargedSpecies::add(size_t n, std::function<void(double&, double&)> sampler) {
    
    m_v.resize(m_n + n);
    m_x.resize(m_n + n);
    m_f.resize(m_n + n);
    
    for(size_t i = m_n; i < m_n + n; i++) {
        sampler(m_v[i], m_x[i]);
        m_f[i] = 0;
    }

    m_n += n;
}