#include <HYPRE_struct_ls.h>
#include <spark/constants/constants.h>

#include "_hypre_utilities.h"
#include "spark/core/matrix.h"
#include "spark/core/vec.h"
#include "spark/electromagnetics/poisson.h"

using namespace spark::electromagnetics;
using namespace spark::core;

namespace {
int stencil_indices[5] = {0, 1, 2, 3, 4};
int stencil_offsets[5][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
constexpr double solver_tolerance = 1e-6;
}  // namespace

StructPoissonSolver::~StructPoissonSolver() = default;

struct StructPoissonSolver::Impl {
    HYPRE_StructGrid hypre_grid_ = nullptr;
    HYPRE_StructStencil hypre_stencil_ = nullptr;
    HYPRE_StructMatrix hypre_A_ = nullptr;
    HYPRE_StructVector hypre_b_ = nullptr;
    HYPRE_StructVector hypre_x_ = nullptr;
    HYPRE_StructSolver hypre_solver_ = nullptr;

    struct BoundaryRef {
        IntVec<2> pos;
        IntVec<2> offset;
        Region* boundary = nullptr;
    };

    struct CellProp {
        CellType cell_type = CellType::Internal;
        Region* region = nullptr;
    };

    TMatrix<CellProp, 2> cells_;
    std::vector<Region> boundaries_;

    DomainProp prop_;

    Impl(const DomainProp& prop, const std::vector<Region>& boundaries);
    void assemble();
    void solve(Matrix<2>& out, const Matrix<2>& rho);
    ~Impl();

private:
    void create_grid();
    void create_matrices();
    void create_stencil();

    void set_cells();
    void set_stencils();

    CellType get_cell(int i, int j);
    core::Matrix<2> input_cache_;
    std::vector<BoundaryRef> boundary_refs_;
};

StructPoissonSolver::Impl::Impl(const DomainProp& prop, const std::vector<Region>& boundaries)
    : boundaries_(boundaries), prop_(prop) {
    input_cache_.resize(prop.extents.to<size_t>());
}

void StructPoissonSolver::Impl::create_grid() {
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &hypre_grid_);

    int lower[] = {0, 0};
    int upper[] = {prop_.extents.x - 1, prop_.extents.y - 1};
    HYPRE_StructGridSetExtents(hypre_grid_, lower, upper);

    HYPRE_StructGridAssemble(hypre_grid_);
}
void StructPoissonSolver::Impl::create_matrices() {
    HYPRE_StructStencilCreate(2, 5, &hypre_stencil_);
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypre_grid_, hypre_stencil_, &hypre_A_);
    HYPRE_StructMatrixInitialize(hypre_A_);

    HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid_, &hypre_b_);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid_, &hypre_x_);
    HYPRE_StructVectorInitialize(hypre_b_);
    HYPRE_StructVectorInitialize(hypre_x_);
}
void StructPoissonSolver::Impl::create_stencil() {
    for (int n = 0; n < 5; n++)
        HYPRE_StructStencilSetElement(hypre_stencil_, n, &stencil_offsets[n][0]);
}

void StructPoissonSolver::Impl::set_cells() {
    cells_.resize(ULongVec<2>(prop_.extents.x, prop_.extents.y));
    cells_.fill({});

    for (const auto& b : boundaries_) {
        cells_.fill({b.region_type, const_cast<Region*>(&b)}, b.lower_left, b.upper_right);
    }
}

void StructPoissonSolver::Impl::set_stencils() {
    const auto [sx, sy] = cells_.size();
    const auto [dx, dy] = prop_.dx;

    const double kx = 1.0 / (dx * dx);
    const double ky = 1.0 / (dy * dy);

    for (int i = 0; i < sx; ++i) {
        for (int j = 0; j < sy; ++j) {
            int index[] = {i, j};

            switch (cells_(i, j).cell_type) {
                case CellType::Internal: {
                    double stencil_internal[] = {-2.0 * (kx + ky), kx, kx, ky, ky};

                    for (int p = 1; p < 5; ++p) {
                        const int pos[] = {i + stencil_offsets[p][0], j + stencil_offsets[p][1]};
                        if (get_cell(pos[0], pos[1]) == CellType::BoundaryDirichlet) {
                            stencil_internal[p] = 0.0;
                            boundary_refs_.push_back(
                                {{i, j},
                                 {stencil_offsets[p][0], stencil_offsets[p][1]},
                                 cells_(pos[0], pos[1]).region});
                        }
                    }

                    HYPRE_StructMatrixSetValues(hypre_A_, index, 5, stencil_indices,
                                                stencil_internal);
                    break;
                }
                case CellType::BoundaryDirichlet: {
                    double stencil_dirichlet[] = {1.0};
                    HYPRE_StructMatrixSetValues(hypre_A_, index, 1, stencil_indices,
                                                stencil_dirichlet);
                    break;
                }
                case CellType::BoundaryNeumann: {
                    double stencil_neumann[] = {-2.0 * (kx + ky), kx, kx, ky, ky};

                    // check x
                    if (get_cell(i - 1, j) == CellType::External) {
                        stencil_neumann[1] = 0.0;
                        stencil_neumann[2] = 2.0 * kx;
                    } else if (get_cell(i + 1, j) == CellType::External) {
                        stencil_neumann[1] = 2.0 * kx;
                        stencil_neumann[2] = 0.0;
                    }

                    // check y
                    if (get_cell(i, j - 1) == CellType::External) {
                        stencil_neumann[3] = 0.0;
                        stencil_neumann[4] = 2.0 * ky;
                    } else if (get_cell(i, j + 1) == CellType::External) {
                        stencil_neumann[3] = 2.0 * ky;
                        stencil_neumann[4] = 0.0;
                    }

                    HYPRE_StructMatrixSetValues(hypre_A_, index, 5, stencil_indices,
                                                stencil_neumann);
                    break;
                }
                default: {
                    break;
                }
            }
        }
    }
}
CellType StructPoissonSolver::Impl::get_cell(int i, int j) {
    if (i >= 0 && i < prop_.extents.x && j >= 0 && j < prop_.extents.y) {
        return cells_(i, j).cell_type;
    }

    return CellType::External;
}

void StructPoissonSolver::Impl::assemble() {
    create_grid();
    set_cells();

    create_matrices();
    create_stencil();
    set_stencils();

    HYPRE_StructMatrixAssemble(hypre_A_);
}
void StructPoissonSolver::Impl::solve(Matrix<2>& out, const Matrix<2>& rho) {
    HYPRE_StructVectorSetConstantValues(hypre_x_, 0.0);
    HYPRE_StructVectorSetConstantValues(hypre_b_, 0.0);

    input_cache_.data().assign(rho.data().begin(), rho.data().end());
    auto& cache = input_cache_.data();
    constexpr double k = -1.0 / constants::eps0;
    for (size_t i = 0; i < input_cache_.count(); ++i) {
        cache[i] *= k;
    }

    {
        int ilower[] = {0, 0};
        int iupper[] = {prop_.extents.x - 1, prop_.extents.y - 1};
        HYPRE_StructVectorSetBoxValues(hypre_b_, ilower, iupper, cache.data());
    }

    const auto [dx, dy] = prop_.dx;

    for (const auto& [region_type, lower_left, upper_right, input] : boundaries_) {
        if (region_type == CellType::BoundaryDirichlet && input) {
            input_cache_.resize({static_cast<size_t>(upper_right.x - lower_left.x + 1),
                                 static_cast<size_t>(upper_right.y - lower_left.y + 1)});
            input_cache_.fill(input());
            int ilower[] = {lower_left.x, lower_left.y};
            int iupper[] = {upper_right.x, upper_right.y};
            HYPRE_StructVectorSetBoxValues(hypre_b_, ilower, iupper, input_cache_.data_ptr());
        }
    }

    const double kx = -1.0 / (dx * dx);
    const double ky = -1.0 / (dy * dy);

    for (auto [pos, offset, boundary] : boundary_refs_) {
        int idx[] = {pos.x, pos.y};
        const double kxy = offset.x != 0 ? kx : ky;
        HYPRE_StructVectorAddToValues(hypre_b_, idx, boundary->input() * kxy);
    }

    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &hypre_solver_);
    HYPRE_StructSMGSetTol(hypre_solver_, solver_tolerance);
    HYPRE_StructSMGSetup(hypre_solver_, hypre_A_, hypre_b_, hypre_x_);
    HYPRE_StructSMGSolve(hypre_solver_, hypre_A_, hypre_b_, hypre_x_);

    {
        out.resize(prop_.extents.to<size_t>());
        for (int i = 0; i < prop_.extents.x; ++i) {
            int idx1[] = {i, 0};
            int idx2[] = {i, prop_.extents.y - 1};
            HYPRE_StructVectorGetBoxValues(hypre_x_, idx1, idx2, &out(i, 0));
        }
    }

    HYPRE_StructSMGDestroy(hypre_solver_);
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
    : impl_(std::make_unique<Impl>(prop, regions)) {
    impl_->assemble();
}

void StructPoissonSolver::solve(core::Matrix<2>& out, const core::Matrix<2>& rho) {
    impl_->solve(out, rho);
}