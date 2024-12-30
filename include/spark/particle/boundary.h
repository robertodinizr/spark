#pragma once

#include <cstdint>

#include "spark/core/matrix.h"
#include "spark/core/vec.h"
#include "spark/particle/species.h"
#include "spark/spatial/grid.h"

namespace spark::particle {

template <unsigned NV>
void apply_symmetric_boundary(ChargedSpecies<1, NV>& species, double xmin, double xmax);

template <unsigned NV>
void apply_absorbing_boundary(ChargedSpecies<1, NV>& species, double xmin, double xmax);

enum class BoundaryType { Specular, Absorbing };

struct TiledBoundary {
    core::IntVec<2> lower_left, upper_right;
    BoundaryType boundary_type;
};

class TiledBoundary2D {
public:
    TiledBoundary2D(const spatial::GridProp<2>& grid_prop,
                    const std::vector<TiledBoundary>& boundaries,
                    double dt);

    void apply(Species<2, 3>* species);

private:
    void add_boundary(const TiledBoundary& boundary, uint8_t id);

    core::TMatrix<uint8_t, 2> cells_;
    spatial::GridProp<2> gprop_;
    std::vector<TiledBoundary> boundaries_;
    double dt_;
};

}  // namespace spark::particle