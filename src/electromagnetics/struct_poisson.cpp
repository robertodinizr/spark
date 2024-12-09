#include <HYPRE_struct_ls.h>

#include "spark/core/matrix.h"
#include "spark/core/vec.h"
#include "spark/electromagnetics/poisson.h"

using namespace spark::electromagnetics;
using namespace spark::core;

namespace {
constexpr int stencil_indices[5] = {0, 1, 2, 3, 4};
constexpr int stencil_offsets[5][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
}  // namespace

StructPoissonSolver::~StructPoissonSolver() = default;

struct StructPoissonSolver::Impl {
    HYPRE_StructGrid hypre_grid_ = nullptr;
    HYPRE_StructStencil hypre_stencil_ = nullptr;
    HYPRE_StructMatrix hypre_A_ = nullptr;
    HYPRE_StructVector hypre_b_ = nullptr;
    HYPRE_StructVector hypre_x_ = nullptr;
    HYPRE_StructSolver hypre_solver_ = nullptr;

    TMatrix<CellType, 2> cells_;
    std::vector<Region> boundaries_;

    DomainProp prop_;

    Impl(const DomainProp& prop, const std::vector<Region>& boundaries);

    void create_grid();
    void create_matrix();

    void set_cells();
    void assemble();
    ~Impl();
};

StructPoissonSolver::Impl::Impl(const DomainProp& prop, const std::vector<Region>& boundaries)
    : boundaries_(boundaries), prop_(prop) {
    assemble();
}

void StructPoissonSolver::Impl::create_grid() {
    HYPRE_StructGridCreate(MPI_COMM_SELF, 2, &hypre_grid_);

    int lower[] = {0, 0};
    int upper[] = {static_cast<int>(prop_.nx) - 1, static_cast<int>(prop_.ny) - 1};
    HYPRE_StructGridSetExtents(hypre_grid_, lower, upper);

    HYPRE_StructGridAssemble(hypre_grid_);
}
void StructPoissonSolver::Impl::create_matrix() {
    HYPRE_StructStencilCreate(2, 5, &hypre_stencil_);
    HYPRE_StructMatrixCreate(MPI_COMM_SELF, hypre_grid_, hypre_stencil_, &hypre_A_);
    HYPRE_StructMatrixInitialize(hypre_A_);
}
void StructPoissonSolver::Impl::set_cells() {
    cells_.resize({prop_.nx, prop_.ny});
    cells_.fill(CellType::Internal);

    const auto size = ULongVec<2>{prop_.nx, prop_.ny};
    for (const auto& [region_type, lower_left, upper_right] : boundaries_) {
        cells_.fill(region_type, lower_left, upper_right);
    }
}

void StructPoissonSolver::Impl::assemble() {
    create_grid();
    set_cells();
    create_matrix();
}

StructPoissonSolver::Impl::~Impl() {
    if (hypre_grid_)
        HYPRE_StructGridDestroy(hypre_grid_);
    if (hypre_stencil_)
        HYPRE_StructStencilDestroy(hypre_stencil_);
    if (hypre_A_)
        HYPRE_StructMatrixDestroy(hypre_A_);
    if (hypre_b_)
        HYPRE_StructVectorDestroy(hypre_b_);
    if (hypre_x_)
        HYPRE_StructVectorDestroy(hypre_x_);
}

StructPoissonSolver::StructPoissonSolver(const StructPoissonSolver::DomainProp& prop,
                                         const std::vector<Region>& regions)
    : impl_(std::make_unique<Impl>(prop, regions)) {}
