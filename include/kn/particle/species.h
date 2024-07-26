#pragma once

#include <vector>
#include <functional>

namespace kn::particle {

    class ChargedSpecies {
    public:
        ChargedSpecies() = default;
        ChargedSpecies(double q, double m) : q(q), m(m) {}
        
        void add(size_t n);
        void add(size_t n, std::function<void(double&, double&)> sampler);
        
        double* v() const { return (double*) vs.data(); };
        double* x() const { return (double*) xs.data(); };
        double* f() const { return (double*) fs.data(); };
        size_t count() const { return n; }

    private:
        std::vector<double> vs;
        std::vector<double> xs;
        std::vector<double> fs;
        size_t n = 0;

        double q = 0;
        double m = 0;
    };
}