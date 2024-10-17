#include "kn/particle/species.h"
#include <cstddef>

void kn::particle::ChargedSpecies1D1V::add(size_t n) {
    
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

void kn::particle::ChargedSpecies1D1V::add(size_t n, std::function<void(double&, double&)> sampler) {
    
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

void kn::particle::ChargedSpecies1D3V::add(size_t n, std::function<void(core::Vec3&, double&)> sampler) {
    
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

template <unsigned NX, unsigned NV>
void kn::particle::ChargedSpecies<NX, NV>::add(size_t n) {
    m_v.resize(m_n + n);
    m_x.resize(m_n + n);
    m_f.resize(m_n + n);

    for(size_t i = m_n; i < m_n + n; i++) {
        m_v[i] = kn::core::Vec<NV>{};
        m_x[i] = kn::core::Vec<NX>{};
        m_f[i] = kn::core::Vec<NX>{};
    }

    m_n += n;
}

template <unsigned NX, unsigned NV>
void kn::particle::ChargedSpecies<NX, NV>::add(size_t n, std::function<void(core::Vec<NV>&, core::Vec<NX>&)> sampler) {
    m_v.resize(m_n + n);
    m_x.resize(m_n + n);
    m_f.resize(m_n + n);
    
    for(size_t i = m_n; i < m_n + n; i++) {
        sampler(m_v[i], m_x[i]);
        m_f[i] = kn::core::Vec<NX>{};
    }

    m_n += n;
}

template <unsigned NX, unsigned NV>
void kn::particle::ChargedSpecies<NX,NV>::add_copy(size_t idx) {
    m_v.push_back(m_v[idx]);
    m_x.push_back(m_x[idx]);
    m_f.push_back(m_f[idx]);
    m_n++;
}

template <unsigned NX, unsigned NV>
void kn::particle::ChargedSpecies<NX,NV>::remove(size_t idx) {
    m_v[idx] = m_v.back();
    m_v.pop_back();

    m_x[idx] = m_x.back();
    m_x.pop_back();

    m_f.pop_back();

    m_n--;
}

template class kn::particle::ChargedSpecies<1,1>;
template class kn::particle::ChargedSpecies<1,3>;
template class kn::particle::ChargedSpecies<2,3>;
template class kn::particle::ChargedSpecies<3,3>;