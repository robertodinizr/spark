#pragma once

#include <memory>
#include <tuple>

#include "spark/constants/constants.h"
#include "spark/random/random.h"
#include "species.h"

namespace spark::particle {

namespace distributions {
inline auto maxwell(const double temperature, const double mass, const double drift = 0) {
    const double vth = std::sqrt(spark::constants::e * temperature / mass);
    return [vth, drift]() { return spark::random::normal(drift, vth); };
}

inline auto uniform(const double vmin, const double vmax) {
    return [vmin, vmax]() { return spark::random::uniform(vmin, vmax); };
}

inline auto single_value(const double val) {
    return [val]() { return val; };
}

inline double maxwellian_flux(const double std, const double drift) {
    return std * sqrt(-2 * log(spark::random::uniform())) + drift;
}

}  // namespace distributions

template <unsigned NX, unsigned NV>
class Emitter {
public:
    Emitter(double rate) : rate_(rate) {}

    using Ref = std::unique_ptr<Emitter<NX, NV>>;
    virtual void emit(Species<NX, NV>& species, double n) = 0;
    void emit(Species<NX, NV>& species) { emit(species, rate_); }
    void set_rate(double rate) { rate_ = rate; }
    virtual ~Emitter() {}

protected:
    double rate_ = 0;
};

template <unsigned NX, unsigned NV, typename... Args>
    requires(sizeof...(Args) == NX + NV)
class IndividualComponentEmitter : public Emitter<NX, NV> {
    using super = Emitter<NX, NV>;

public:
    IndividualComponentEmitter(double rate, Args&&... distributions)
        : super(rate), distributions_(std::make_tuple(std::forward<Args>(distributions)...)) {}

    void emit(Species<NX, NV>& species, double n) override {
        double fn = floor(n);
        int np = spark::random::uniform() <= (n - fn) ? fn + 1 : fn;

        species.add(np, [this](auto& v, auto& x) {
            if constexpr (NX > 0)
                x.x = std::get<0>(distributions_)();
            if constexpr (NX > 1)
                x.y = std::get<1>(distributions_)();
            if constexpr (NX > 2)
                x.z = std::get<2>(distributions_)();

            if constexpr (NV > 0)
                v.x = std::get<NX + 0>(distributions_)();
            if constexpr (NV > 1)
                v.x = std::get<NX + 1>(distributions_)();
            if constexpr (NV > 2)
                v.x = std::get<NX + 2>(distributions_)();
        });
    }

private:
    std::tuple<Args...> distributions_;
};

template <unsigned NX, unsigned NV, typename F>
class CombinedComponentEmitter : public Emitter<NX, NV> {
    using super = Emitter<NX, NV>;

public:
    CombinedComponentEmitter(double rate, F&& distribution)
        : super(rate), distribution_(std::forward<F>(distribution)) {}

    void emit(Species<NX, NV>& species, double n) override {
        double fn = floor(n);
        int np = spark::random::uniform() <= (n - fn) ? fn + 1 : fn;

        species.add(np, [this](auto& v, auto& x) {
            const auto sample = distribution_();
            x = std::get<0>(sample);
            v = std::get<1>(sample);
        });
    }

private:
    F distribution_;
};

template <unsigned NX, unsigned NV, typename... Args>
    requires(sizeof...(Args) == NX + NV)
inline Emitter<NX, NV>::Ref make_emitter(Args&&... individual_distributions) {
    return std::make_unique<IndividualComponentEmitter<NX, NV, Args...>>(
        0, std::forward<Args>(individual_distributions)...);
}

template <unsigned NX, unsigned NV, typename... Args>
    requires(sizeof...(Args) == NX + NV)
inline Emitter<NX, NV>::Ref make_emitter(double rate, Args&&... individual_distributions) {
    return std::make_unique<IndividualComponentEmitter<NX, NV, Args...>>(
        rate, std::forward<Args>(individual_distributions)...);
}

template <unsigned NX, unsigned NV, typename F>
inline Emitter<NX, NV>::Ref make_emitter(F&& combined_distribution) {
    return std::make_unique<CombinedComponentEmitter<NX, NV, F>>(
        0, std::forward<F>(combined_distribution));
}

template <unsigned NX, unsigned NV, typename F>
inline Emitter<NX, NV>::Ref make_emitter(double rate, F&& combined_distribution) {
    return std::make_unique<CombinedComponentEmitter<NX, NV, F>>(
        rate, std::forward<F>(combined_distribution));
}

}  // namespace spark::particle