#include "kn/particle/species.h"
#include <cstddef>

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

void kn::particle::ChargedSpecies1D3V::add(size_t n) {
    
    m_v.resize(m_n + n);
    m_x.resize(m_n + n);
    m_f.resize(m_n + n);
    
    for(size_t i = m_n; i < m_n + n; i++) {
        m_v[i] = {};
        m_x[i] = 0.0;
        m_f[i] = {};
    }

    m_n += n;
}

void kn::particle::ChargedSpecies1D3V::add(size_t n, std::function<void(ChargedSpecies1D3V::Vec3&, double&)> sampler) {
    
    m_v.resize(m_n + n);
    m_x.resize(m_n + n);
    m_f.resize(m_n + n);
    
    for(size_t i = m_n; i < m_n + n; i++) {
        sampler(m_v[i], m_x[i]);
        m_f[i] = 0.0;
    }

    m_n += n;
}

void kn::particle::ChargedSpecies1D3V::add_copy(size_t idx) {
    m_v.push_back(m_v[idx]);
    m_x.push_back(m_x[idx]);
    m_f.push_back(m_f[idx]);
    m_n++;
}

void kn::particle::ChargedSpecies1D3V::remove(size_t idx) {
    m_v[idx] = m_v.back();
    m_v.pop_back();

    m_x[idx] = m_x.back();
    m_x.pop_back();

    m_f.pop_back();

    m_n--;
}