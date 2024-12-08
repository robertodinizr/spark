#include <HYPRE_struct_ls.h>

#include "spark/electromagnetics/poisson.h"

using namespace spark::electromagnetics;

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

    DomainProp prop_;

    void assemble_grid() {
        HYPRE_StructGridCreate(MPI_COMM_SELF, 2, &hypre_grid_);
        HYPRE_StructGridSetExtents(hypre_grid_, inner_ll.val, inner_ur.val);
        n_solve +=
            (inner_ur.val[0] - inner_ll.val[0] + 1) * (inner_ur.val[1] - inner_ll.val[1] + 1);

        // Neumann boundary points
        for (int i = 0; i < n_neumann; i++) {
            ll = {neumann_boxes.val[i * 4 + 0], neumann_boxes.val[i * 4 + 1]};
            ur = {neumann_boxes.val[i * 4 + 2], neumann_boxes.val[i * 4 + 3]};

            HYPRE_StructGridSetExtents(hypre_grid_, ll.val, ur.val);
            n_solve += (ur.val[0] - ll.val[0] + 1) * (ur.val[1] - ll.val[1] + 1);
        }
        HYPRE_StructGridAssemble(hypre_grid_);
    }

    void assemble(const std::vector<Boundary>& boundaries) { assemble_grid(); }

    ~Impl() {
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
};

StructPoissonSolver::StructPoissonSolver(const StructPoissonSolver::DomainProp& prop,
                                         const std::vector<Boundary>& boundaries)
    : impl_(new Impl{.prop_ = prop}) {
    impl_->assemble(boundaries);
}
