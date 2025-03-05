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
}  // namespace distributions

template <unsigned NX, unsigned NV>
class Emitter {
public:
    using Ref = std::unique_ptr<Emitter<NX, NV>>;
    virtual void emit(Species<NX, NV>& species, int n) = 0;
    void emit(Species<NX, NV>& species) { emit(species, rate_); }
    void set_rate(int rate) { rate_ = rate; }
    virtual ~Emitter() {}

protected:
    int rate_ = 0;
};

template <unsigned NX, unsigned NV, typename... Args>
    requires(sizeof...(Args) == NX + NV)
class EmitterImpl : public Emitter<NX, NV> {
    using super = Emitter<NX, NV>;

public:
    EmitterImpl() = default;
    EmitterImpl(int rate, Args&&... distributions)
        : distributions_(std::make_tuple(std::forward<Args>(distributions)...)) {
        Emitter<NX, NV>::rate_ = rate;
    }

    void emit(Species<NX, NV>& species, int n) override {
        species.add(n, [this](auto& v, auto& x) {
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

template <unsigned NX, unsigned NV, typename... Args>
inline Emitter<NX, NV>::Ref make_emitter(Args&&... distributions) {
    return std::make_unique<EmitterImpl<NX, NV, Args...>>(0, std::forward<Args>(distributions)...);
}

template <unsigned NX, unsigned NV, typename... Args>
inline Emitter<NX, NV>::Ref make_emitter(int rate, Args&&... distributions) {
    return std::make_unique<EmitterImpl<NX, NV, Args...>>(rate,
                                                          std::forward<Args>(distributions)...);
}

}  // namespace spark::particle