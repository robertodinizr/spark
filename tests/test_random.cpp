#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <vector>
#include <numeric>
#include <cmath>

#include "spark/random/random.h"

using namespace spark::random;

TEST_CASE("Random number generation") {

    const uint64_t seed = 123456789;

    SECTION("Initialize") {
        std::vector<double> sequence1, sequence2;
        const int sequence_size = 100;
        sequence1.reserve(sequence_size);
        sequence2.reserve(sequence_size);

        initialize(seed);
        for (int i = 0; i < sequence_size; ++i) {
            sequence1.push_back(uniform());
        }

        initialize(seed);
        for (int i = 0; i < sequence_size; ++i) {
            sequence2.push_back(uniform());
        }

        REQUIRE(sequence1.size() == sequence_size);
        REQUIRE(sequence2.size() == sequence_size);

        REQUIRE(sequence1 == sequence2);
    }

    SECTION("Uniform distribution") {
        initialize(seed);
        const int num_samples = 10000;

        for (int i = 0; i < num_samples; ++i) {
            double val = uniform();
            REQUIRE(val >= 0.0);
            REQUIRE(val <= 1.0);
        }

        const double vmax = 50.0;
        for (int i = 0; i < num_samples; ++i) {
            double val = uniform(vmax);
            REQUIRE(val >= 0.0);
            REQUIRE(val <= vmax);
        }

        const double vmin = -25.0;
        for (int i = 0; i < num_samples; ++i) {
            double val = uniform(vmin, vmax);
            REQUIRE(val >= vmin);
            REQUIRE(val <= vmax);
        }
    }

    SECTION("Normal distribution") {
        initialize(seed);
        const int num_samples = 20000;
        const double mean = 10.0;
        const double std_dev = 2.5;

        std::vector<double> samples;
        samples.reserve(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            samples.push_back(normal(mean, std_dev));
        }

        double sample_mean = std::accumulate(samples.begin(), samples.end(), 0.0) / num_samples;

        double sum_sq_diff = 0.0;
        for(const auto& s : samples) {
            sum_sq_diff += (s - sample_mean) * (s - sample_mean);
        }
        double sample_std_dev = std::sqrt(sum_sq_diff / (num_samples - 1));

        REQUIRE(sample_mean == Catch::Approx(mean).margin(0.1));
        REQUIRE(sample_std_dev == Catch::Approx(std_dev).margin(0.1));
    }
}