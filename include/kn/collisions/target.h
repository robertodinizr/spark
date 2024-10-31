#pragma once

#include "kn/core/vec.h"

namespace kn::collisions {

template <unsigned NX, unsigned NV>
class Target {
public:
    virtual double dens_at(const core::Vec<NX>& pos) = 0;
    virtual double dens_max() = 0;
    virtual double temperature() = 0;
    virtual ~Target(){};
};

template <unsigned NX, unsigned NV>
class StaticUniformTarget : public Target<NX, NV> {
public:
    StaticUniformTarget<NX, NV>(double density, double temperature)
        : m_density(density), m_temperature(temperature) {}

    virtual double dens_at(const core::Vec<NX>& pos) override { return m_density; }
    virtual double dens_max() override { return m_density; }
    virtual double temperature() override { return m_temperature; }

private:
    double m_density = 0.0;
    double m_temperature = 0.0;
};
}  // namespace kn::collisions
