#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "spark/particle/emitter.h"
#include "spark/particle/species.h"
#include "spark/random/random.h"
#include <numeric>

using namespace spark::particle;
using namespace spark::core;

TEST_CASE("Emitter Distributions") {
    spark::random::initialize(12345);

    SECTION("single_value distribution") {
        auto dist = distributions::single_value(10.5);
        REQUIRE(dist() == Catch::Approx(10.5));
        REQUIRE(dist() == Catch::Approx(10.5));
    }

    SECTION("uniform distribution") {
        auto dist = distributions::uniform(-5.0, 5.0);
        for (int i=0; i<100; ++i) {
            double val = dist();
            REQUIRE(val >= -5.0);
            REQUIRE(val <= 5.0);
        }
    }
}

TEST_CASE("IndividualComponentEmitter Tests") {
    spark::random::initialize(12345);
    ChargedSpecies<1, 3> species(1.0, 1.0);

    auto emitter = make_emitter<1, 3>(
        distributions::single_value(5.0),
        distributions::single_value(100.0),
        distributions::single_value(200.0),
        distributions::single_value(300.0)
    );

    emitter->emit(species, 2.0);

    REQUIRE(species.n() == 2);

    for (size_t i = 0; i < species.n(); ++i) {
        REQUIRE(species.x()[i].x == Catch::Approx(5.0));
        REQUIRE(species.v()[i].x == Catch::Approx(100.0));
        REQUIRE(species.v()[i].y == Catch::Approx(200.0));
        REQUIRE(species.v()[i].z == Catch::Approx(300.0));
    }
}


TEST_CASE("CombinedComponentEmitter Tests") {
    spark::random::initialize(12345);
    ChargedSpecies<2, 3> species(1.0, 1.0);

    auto combined_dist = []() {
        return std::make_pair(Vec<2>{1.0, 2.0}, Vec<3>{10.0, 20.0, 30.0});
    };
    auto emitter = make_emitter<2, 3>(1.0, combined_dist);

    emitter->emit(species, 1.0);

    REQUIRE(species.n() == 1);
    REQUIRE(species.x()[0].x == Catch::Approx(1.0));
    REQUIRE(species.x()[0].y == Catch::Approx(2.0));
    REQUIRE(species.v()[0].x == Catch::Approx(10.0));
    REQUIRE(species.v()[0].y == Catch::Approx(20.0));
    REQUIRE(species.v()[0].z == Catch::Approx(30.0));
}

TEST_CASE("Emitter Fractional Rate") {
    spark::random::initialize(12345);
    ChargedSpecies<1, 1> species(1.0, 1.0);
    auto emitter = make_emitter<1, 1>(distributions::single_value(0.0), distributions::single_value(0.0));

    int n_runs = 1000;
    double rate = 0.7;
    int particles_created = 0;

    for (int i = 0; i < n_runs; ++i) {
        size_t n_before = species.n();
        emitter->emit(species, rate);
        size_t n_after = species.n();
        particles_created += (n_after - n_before);
    }
    
    REQUIRE(static_cast<double>(particles_created) / n_runs == Catch::Approx(rate).margin(0.05));
}