#pragma once

#include <vector>

#include "spark/core/vec.h"

namespace spark::particle {

template <unsigned NX, unsigned NV>
class Species {
public:
    Species() = default;
    explicit Species(const double m) : m_(m) {}

    double m() const { return m_; }
    size_t n() const { return x_.size(); }
    core::Vec<NV>* v() const { return const_cast<core::Vec<NV>*>(v_.data()); };
    core::Vec<NX>* x() const { return const_cast<core::Vec<NX>*>(x_.data()); };

    void add(size_t n) {
        const auto n_current = x_.size();
        v_.resize(n_current + n);
        x_.resize(n_current + n);
    }

    void add(size_t n, auto sampler) {
        const auto n_current = x_.size();

        v_.resize(n_current + n);
        x_.resize(n_current + n);

        for (size_t i = n_current; i < n_current + n; ++i) {
            sampler(v_[i], x_[i]);
        }
    }

    template <auto SamplerFunc(core::Vec<NV>&, core::Vec<NX>&)->void>
    void add(size_t n) {
        const auto n_current = x_.size();

        v_.resize(n_current + n);
        x_.resize(n_current + n);

        for (size_t i = n_current; i < n_current + n; ++i) {
            SamplerFunc(v_[i], x_[i]);
        }
    }

    void add_copy(size_t idx) {
        v_.push_back(v_[idx]);
        x_.push_back(x_[idx]);
    }

    void remove(size_t idx) {
        v_[idx] = v_.back();
        v_.pop_back();

        x_[idx] = x_.back();
        x_.pop_back();
    }

private:
    std::vector<core::Vec<NV>> v_;
    std::vector<core::Vec<NX>> x_;
    double m_ = 0;
};

template <unsigned NX, unsigned NV>
class ChargedSpecies : public Species<NX, NV> {
public:
    ChargedSpecies() = default;
    ChargedSpecies(double q, double m) : Species<NX, NV>(m), q_(q) {}
    double q() const { return q_; }

private:
    double q_ = 0;
};

}  // namespace spark::particle