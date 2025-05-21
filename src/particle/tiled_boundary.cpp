
#include <parallel_hashmap/phmap.h>

#include <cstddef>
#include <cstdint>
#include <queue>

#include "spark/core/vec.h"
#include "spark/particle/boundary.h"
#include "spark/spatial/grid.h"

#define SPARK_TILED_BOUNDARY_CHECK_IF_INSIDE

using namespace spark::core;

template <>
struct std::hash<IntVec<2>> {
    std::size_t operator()(const IntVec<2>& k) const {
        using std::hash;
        using std::size_t;
        const auto h1 = hash<int>()(k.x);
        return h1 ^ (hash<int>()(k.y) + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

namespace {

struct CollisionHit {
    spark::core::Vec<2> normal;
    spark::core::Vec<2> pos;
    uint8_t val = 0;
};

#define CMOD(idx, size) ((size + (idx % size)) % size)
#define CMODV(idx, size_vec) CMOD(idx.x, size_vec.x), CMOD(idx.y, size_vec.y)

constexpr size_t padding_ = 4;

int cmod(int idx, int size) {
    return ((size + (idx % size)) % size);
}

#define FLOOR_F(x, fp) \
    int x = int(fp);   \
    if (x > fp)        \
        x--;

// Assuming cell size is 1x1
void grid_raycast(const spark::core::TMatrix<uint8_t, 2>& grid,
                  const Vec<2>& a,
                  const Vec<2>& b,
                  CollisionHit& hit) {
    hit.val = 0;
    FLOOR_F(current_index_x, a.x)
    FLOOR_F(current_index_y, a.y)
    FLOOR_F(end_index_x, b.x)
    FLOOR_F(end_index_y, b.y)

    const auto sz = grid.size().to<int>();
    auto dir = (b - a).normalized();

    Vec<2> t{0.0, 0.0};
    IntVec<2> step = {0, 0};
    Vec<2> k = {std::abs(1.0 / dir.x), std::abs(1.0 / dir.y)};

    if (dir.x < 0) {
        step.x = -1;
        t.x = (a.x - static_cast<double>(current_index_x)) * k.x;
    } else {
        step.x = 1;
        t.x = (static_cast<double>(current_index_x + 1) - a.x) * k.x;
    }

    if (dir.y < 0) {
        step.y = -1;
        t.y = (a.y - static_cast<double>(current_index_y)) * k.y;
    } else {
        step.y = 1;
        t.y = (static_cast<double>(current_index_y + 1) - a.y) * k.y;
    }

    double distance = 0.0;
    Vec<2> normal = {0.0, 0.0};

    while (current_index_x != end_index_x || current_index_y != end_index_y) {
        if (t.x < t.y) {
            current_index_x += step.x;
            distance = t.x;
            t.x += k.x;
            normal = {-(double)step.x, 0};
        } else {
            current_index_y += step.y;
            distance = t.y;
            t.y += k.y;
            normal = {0, -(double)step.y};
        }

        const int ki = CMOD(current_index_x, sz.x);
        const int kj = CMOD(current_index_y, sz.y);

        if (uint8_t val = grid(ki, kj)) {
            hit.normal = normal;
            hit.pos = a + dir * distance;
            hit.val = val;
            return;
        }
    }
}

inline Vec<2> reflect(const Vec<2>& x, const Vec<2>& n, const Vec<2>& contact) {
    // Omitting terms:
    // 1) 2 * n.x * n.y * (contact.y - x.y) in x
    // 2) 2 * n.x * n.y * (contact.x - x.x) in y
    // Since they will always be zero if n is axis aligned
    return {2.0 * (contact.x - x.x) * n.x * n.x + x.x, 2.0 * (contact.y - x.y) * n.y * n.y + x.y};
}

inline Vec<3> reflect(const Vec<3>& v, const Vec<2>& n) {
    return {v.x * (1.0 - 2.0 * n.x * n.x), v.y * (1.0 - 2.0 * n.y * n.y), v.z};
}
}  // namespace

namespace spark::particle {
TiledBoundary2D::TiledBoundary2D(const spatial::GridProp<2>& grid_prop,
                                 const std::vector<TiledBoundary>& boundaries,
                                 const double dt,
                                 const bool empty_box)
    : gprop_(grid_prop), boundaries_(boundaries), dt_(dt), empty_box_(empty_box) {
    cells_.resize(grid_prop.n + padding_);
    distance_cells_.resize(grid_prop.n);

    sx_ = gprop_.n.x - 1;
    sy_ = gprop_.n.y - 1;
    for (uint8_t i = 0; i < boundaries_.size(); ++i) {
        // Add circular_mod
        add_boundary(boundaries_[i], i + 1);
    }

    set_distance_cells();
}

uint8_t TiledBoundary2D::cell(int i, int j) const {
    const auto sz = cells_.size().to<int>();
    return cells_(CMOD(i, sz.x), CMOD(j, sz.y));
}

uint8_t TiledBoundary2D::cell(const Vec<2>& pos) const {
    auto idx = (pos / gprop_.dx).apply<std::floor>().to<int>();
    const auto sz = cells_.size().to<int>();
    return cells_(CMOD(idx.x, sz.x), CMOD(idx.y, sz.y));
}

void TiledBoundary2D::add_boundary(const TiledBoundary& b, uint8_t id) {
    const auto imin = std::min(b.lower_left.x, b.upper_right.x);
    const auto imax = std::max(b.lower_left.x, b.upper_right.x);
    const auto jmin = std::min(b.lower_left.y, b.upper_right.y);
    const auto jmax = std::max(b.lower_left.y, b.upper_right.y);
    const auto sz = cells_.size().to<int>();

    for (int i = imin; i <= imax; ++i)
        for (int j = jmin; j <= jmax; ++j) {
            int ki = CMOD(i, sz.x);
            int kj = CMOD(j, sz.y);

            cells_(ki, kj) = id;
        }
}

int TiledBoundary2D::bfs_closest_boundary(int i, int j) {
    phmap::flat_hash_set<IntVec<2>> visited;
    std::queue<IntVec<2>> to_visit;

    to_visit.push(IntVec<2>{i, j});
    while (!to_visit.empty()) {
        auto current = to_visit.front();
        to_visit.pop();

        if (cell(current.x, current.y)) {
            return std::abs(current.x - i) + std::abs(current.y - j);
        }

        if (auto n = IntVec<2>{current.x + 1, current.y}; !visited.contains(n)) {
            to_visit.push(n);
            visited.insert(n);
        }
        if (auto n = IntVec<2>{current.x - 1, current.y}; !visited.contains(n)) {
            to_visit.push(n);
            visited.insert(n);
        }
        if (auto n = IntVec<2>{current.x, current.y + 1}; !visited.contains(n)) {
            to_visit.push(n);
            visited.insert(n);
        }
        if (auto n = IntVec<2>{current.x, current.y - 1}; !visited.contains(n)) {
            to_visit.push(n);
            visited.insert(n);
        }
    }

    return 0;
}

void TiledBoundary2D::set_distance_cells() {
    for (int i = 0; i < sx_; i++) {
        for (int j = 0; j < sy_; j++) {
            distance_cells_(i, j) = bfs_closest_boundary(i, j);
        }
    }
}

void TiledBoundary2D::apply(Species<2, 3>& species, Stats* stats) {
    int n = species.n();
    auto* x = species.x();
    auto* v = species.v();

    if (stats) {
        stats->ncolls.resize(boundaries_.size());
        std::fill(stats->ncolls.begin(), stats->ncolls.end(), 0);
    }

    for (int i = 0; i < n; ++i) {
        auto& x1 = x[i];

        if (empty_box_ && (x1.x > 0 && x1.x < gprop_.l.x && x1.y > 0 && x1.y < gprop_.l.y))
            continue;

        auto& v1 = v[i];
        const auto x0 = core::Vec<2>{x1.x - v1.x * dt_, x1.y - v1.y * dt_};

        auto x0_tmp = x0 / gprop_.dx;
        auto x1_tmp = x1 / gprop_.dx;

        const auto x0_cell = x0_tmp.apply<std::floor>().to<int>();
        const auto x1_cell = x1_tmp.apply<std::floor>().to<int>();

        const uint8_t distance_to_boundary = distance_cells_(x0_cell.x, x0_cell.y);
        if (std::abs(x0_cell.x - x1_cell.x) + std::abs(x0_cell.y - x1_cell.y) <
            distance_to_boundary)
            continue;

#ifdef SPARK_TILED_BOUNDARY_CHECK_IF_INSIDE
        // If particle is inside wall, remove it to avoid errors
        if (distance_to_boundary == 0) {
            species.remove(i);
            i--;
            n--;
            continue;
        }
#endif

        CollisionHit hit{0};

        while (true) {
            grid_raycast(cells_, x0_tmp, x1_tmp, hit);

            // val == 0 means that no boundary was found
            if (!hit.val)
                break;

            const auto& b = boundaries_[hit.val - 1];
            const auto btype = b.boundary_type;

            if (stats) {
                stats->ncolls[hit.val - 1]++;
            }

            if (btype == BoundaryType::Absorbing) {
                // TODO(lui): Check if this is OK
                species.remove(i);
                i--;  // check ith particle again since the particle is replaced during removal
                n--;  // decrease the number of particles
                break;
            } else if (btype == BoundaryType::Specular) {
                // Specular reflection
                x1 = reflect(x1, hit.normal, hit.pos * gprop_.dx);
                v1 = reflect(v1, hit.normal);

                x1_tmp = x1 / gprop_.dx;
                x0_tmp = hit.pos;
            }

            // TODO(lui): Implement diffuse (Lambertian) reflection
        }
    }
}
}  // namespace spark::particle
