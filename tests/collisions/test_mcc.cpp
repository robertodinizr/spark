#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "spark/collisions/mcc.h"
#include "spark/collisions/basic_reactions.h"
#include "spark/particle/species.h"
#include "spark/random/random.h"
#include <memory>

using namespace spark::collisions;
using namespace spark::particle;
using namespace spark::core;


TEST_CASE("MCCReactionSet Core Logic") {
    spark::random::initialize(12345);

    ChargedSpecies<2, 3> projectile(-spark::constants::e, spark::constants::m_e);
    
    ReactionConfig<2, 3> config;
    config.dt = 1e-7;
    config.dyn = RelativeDynamics::FastProjectile;

    config.reactions.push_back(std::make_unique<reactions::ChargeExchangeCollision<2,3>>(
        reactions::BasicCollisionConfig{}, 
        CrossSection{0.0, {0.1, 1000.0}, {1e-19, 1e-19}}
    ));

    SECTION("No collisions occur when target density is zero") {
        
        config.target = std::make_shared<StaticUniformTarget<2, 3>>(0.0, 300.0);

        projectile.add(1, [](auto& v, auto& x){ v = {1e5, 1e5, 1e5}; });
        
        auto initial_velocity = projectile.v()[0];

        MCCReactionSet<2, 3> mcc_set(&projectile, std::move(config));
        mcc_set.react_all();

        auto final_velocity = projectile.v()[0];
        REQUIRE(final_velocity.x == initial_velocity.x);
        REQUIRE(final_velocity.y == initial_velocity.y);
        REQUIRE(final_velocity.z == initial_velocity.z);
        REQUIRE(projectile.n() == 1);
    }
}