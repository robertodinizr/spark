#include <HYPRE_struct_ls.h>
#include <spark/constants/constants.h>
#include "HYPRE_utilities.h"
#include "_hypre_utilities.h"
#include "log/log.h"
#include "spark/core/matrix.h"
#include "spark/core/vec.h"
#include "spark/em/poisson.h"

using namespace spark::core;
using namespace spark::em;
using namespace spark::spatial;

namespace {
    int stencil_indices[5] = {0, 1, 2, 3, 4};
    int opposite_indices[] = {0, 2, 1, 4, 3};
    int stencil_offsets[5][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    constexpr double solver_tolerance = 1e-6;
    bool hypre_initialized = false;
}

CylindricalPoissonSolver2D::~CylindricalPoissonSolver2D() = default;

struct CylindricalPoissonSolver2D::Impl {
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

    void set_stencils();
    void set_cells();

    CellType get_cell(int i, int j);
    core::Matrix<2> input_cache_;
    std::vector<BoundaryRef> boundary_refs_;
};

CylindricalPoissonSolver2D::Impl::Impl(const DomainProp& prop, const std::vector<Region>& boundaries)
    : boundaries_(boundaries), prop_(prop) {
    if (!hypre_initialized) {
        HYPRE_Initialize();
        hypre_initialized = true;
    }
    input_cache_.resize(prop.extents.to<size_t>());
}

void CylindricalPoissonSolver2D::Impl::create_grid() {
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &hypre_grid_);
    int lower[] = {0, 0};
    int upper[] = {prop_.extents.x - 1, prop_.extents.y - 1};
    HYPRE_StructGridSetExtents(hypre_grid_, lower, upper);
    HYPRE_StructGridAssemble(hypre_grid_);
}

void CylindricalPoissonSolver2D::Impl::create_matrices() {
    HYPRE_StructStencilCreate(2, 5, &hypre_stencil_);
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypre_grid_, hypre_stencil_, &hypre_A_);
    HYPRE_StructMatrixInitialize(hypre_A_);

    HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid_, &hypre_b_);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid_, &hypre_x_);
    HYPRE_StructVectorInitialize(hypre_b_);
    HYPRE_StructVectorInitialize(hypre_x_);
}

void CylindricalPoissonSolver2D::Impl::create_stencil() {
    for (int n = 0; n < 5; n++)
        HYPRE_StructStencilSetElement(hypre_stencil_, n, &stencil_offsets[n][0]);
}

void CylindricalPoissonSolver2D::Impl::set_cells() {
    cells_.resize(ULongVec<2>(prop_.extents.x, prop_.extents.y));
    cells_.fill({});

    for (const auto& b : boundaries_) {
        cells_.fill({b.region_type, const_cast<Region*>(&b)}, b.lower_left, b.upper_right);
    }
}

void CylindricalPoissonSolver2D::Impl::set_stencils() {
    const auto [sz, sr] = cells_.size();
    const auto [dz, dr] = prop_.dx;
    const double idz = 1.0 / dz;
    const double idr = 1.0 / dr;
    const double idz2 = 1.0 / (dz * dz);
    const double idr2 = 1.0 / (dr * dr);
    const double opposite_k[] = {-(2.0 * idz2 + 2.0 * idr2), idz2, idz2, 0.0,2.0 * idr2};

    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sr; ++j) {
            int index[] = {i, j};

            auto cell = cells_(i, j).cell_type;
            double coeff_center, coeff_right, coeff_left, coeff_up, coeff_down;
            if (cell == CellType::BoundaryDirichlet) {
                double stencil_dirichlet[] = {1.0};
                HYPRE_StructMatrixSetValues(hypre_A_, index, 1, stencil_indices, stencil_dirichlet);
                continue;
            }
            coeff_down = idr2 - 0.5 * idr;
            coeff_left = idz2;
            coeff_center = -2.0 * (idz2 + idr2);
            coeff_right = idz2;
            coeff_up = idr2 + 0.5 * idr;

            double stencil[] = {coeff_center, coeff_left, coeff_right, coeff_down, coeff_up};

            for (int p = 1; p < 5; ++p) {
                const int neighbor_pos[] = {i + stencil_offsets[p][0], j + stencil_offsets[p][1]};
                const auto neighbor_type = get_cell(neighbor_pos[0], neighbor_pos[1]);

                if (neighbor_type == CellType::BoundaryDirichlet) {
                    stencil[p] = 0.0;
                    boundary_refs_.push_back({{i, j},
                    {stencil_offsets[p][0], stencil_offsets[p][1]},
                    cells_(neighbor_pos[0], neighbor_pos[1]).region});
                }

                if (cell == CellType::BoundaryNeumann && neighbor_type == CellType::External) {
                    stencil[p] = 0.0;
                    stencil[stencil_indices[p]] = 2.0 * opposite_k[p];
                }

                if (cell == CellType::Internal && neighbor_type == CellType::External) {
                    SPARK_LOG_WARN("internal node [%d, %d] at boundary with offset [%d, %d]!", i, j,
                                   stencil_offsets[p][0], stencil_offsets[p][1]);
                }
            }
            HYPRE_StructMatrixSetValues(hypre_A_, index, 5, stencil_indices, stencil);
        }
    }
}

CellType CylindricalPoissonSolver2D::Impl::get_cell(int i, int j) {
    if (i >= 0 && i < prop_.extents.x && j >= 0 && j < prop_.extents.y) {
        return cells_(i, j).cell_type;
    }
    return CellType::External;
}

void CylindricalPoissonSolver2D::Impl::assemble() {
    create_grid();
    set_cells();
    create_matrices();
    create_stencil();
    set_stencils();
    HYPRE_StructMatrixAssemble(hypre_A_);
}

void CylindricalPoissonSolver2D::Impl::solve(core::Matrix<2>& out, const core::Matrix<2>& rho) {
    constexpr double k = -1.0 / constants::eps0;
   
    for (int i = 0; i < rho.size().x; i++) {
        for (int j = 0; j < rho.size().y; ++j) {
            int pos[] = {i, j};
            HYPRE_StructVectorSetValues(hypre_b_, pos, rho(i, j) * k);    
	}
    }

    const auto [dz, dr] = prop_.dx;
    for (const auto& [region_type, lower_left, upper_right, input] : boundaries_) {
	if (region_type == CellType::BoundaryDirichlet && input) {
            input_cache_.resize({static_cast<size_t>(upper_right.x - lower_left.x + 1),                                 static_cast<size_t>(upper_right.y - lower_left.y + 1)});
	    input_cache_.fill(input());
	    int ilower[] = {lower_left.x, lower_left.y};
	    int iupper[] = {upper_right.x, upper_right.y};
	    HYPRE_StructVectorSetBoxValues(hypre_b_, ilower, iupper, input_cache_.data_ptr());
	}
    }

    const double kz = -1.0 / (dz * dz);
    const double kr = -1.0 / (dr * dr);

    for (auto [pos, offset, boundary] : boundary_refs_) {
	int idx[] = {pos.x, pos.y};
        const double kxy = offset.x != 0 ? kz : kr;
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

CylindricalPoissonSolver2D::Impl::~Impl() {
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
CylindricalPoissonSolver2D::CylindricalPoissonSolver2D() : impl_(nullptr) {}
CylindricalPoissonSolver2D::CylindricalPoissonSolver2D(const CylindricalPoissonSolver2D::DomainProp& prop,const std::vector<Region>& regions)
: impl_(std::make_unique<Impl>(prop, regions)) {
	impl_->assemble();
}
CylindricalPoissonSolver2D& CylindricalPoissonSolver2D::operator=(CylindricalPoissonSolver2D&& other) noexcept {
    impl_ = std::move(other.impl_);
    return *this;
}
CylindricalPoissonSolver2D::CylindricalPoissonSolver2D(CylindricalPoissonSolver2D&& other) noexcept
: impl_(std::move(other.impl_)) {}
void CylindricalPoissonSolver2D::solve(core::Matrix<2>& out, const core::Matrix<2>& rho) {
    impl_->solve(out, rho);
}
