/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/DebuggingUtilities.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/ibtk_utilities.h"

#include "LS/LSCutCellLaplaceOperator.h"
#include "LS/utility_functions.h"

#include "CellDataFactory.h"
#include "CellVariable.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <ostream>
#include <string>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace LS
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 2;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "CONSERVATIVE_LINEAR_REFINE";
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

static Timer* t_compute_helmholtz;
static Timer* t_extrapolate;
static Timer* t_apply;
static Timer* t_find_cell_centroid;
static Timer* t_find_system;
static Timer* t_solve_system;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

LSCutCellLaplaceOperator::LSCutCellLaplaceOperator(const std::string& object_name,
                                                   Pointer<Database> input_db,
                                                   const bool homogeneous_bc)
    : LaplaceOperator(object_name, homogeneous_bc), d_Q_var(new CellVariable<NDIM, double>(d_object_name + "::ScrVar"))
{
    // Setup the operator to use default scalar-valued boundary conditions.
    setPhysicalBcCoef(nullptr);

    d_robin_bdry = input_db->getBoolWithDefault("robin_boundary", d_robin_bdry);

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_Q_scr_idx = var_db->registerVariableAndContext(
        d_Q_var, var_db->getContext(d_object_name + "::SCRATCH"), input_db->getInteger("stencil_size"));

    LS_DO_ONCE(t_compute_helmholtz =
                   TimerManager::getManager()->getTimer("LS::LSCutCellLaplaceOperator::computeHelmholtzAction()");
               t_extrapolate =
                   TimerManager::getManager()->getTimer("LS::LSCutCellLaplaceOperator::extrapolateToCellCenters()");
               t_apply = TimerManager::getManager()->getTimer("LS::LSCutCellLaplaceOperator::apply()");
               t_find_cell_centroid =
                   TimerManager::getManager()->getTimer("LS::LSCutCellLaplaceOperator::find_cell_centroid");
               t_find_system = TimerManager::getManager()->getTimer("LS::LSCutCellLaplaceOperator::find_system");
               t_solve_system = TimerManager::getManager()->getTimer("LS::LSCutCellLaplaceOperator::solve_system"););
    return;
} // LSCutCellLaplaceOperator()

LSCutCellLaplaceOperator::~LSCutCellLaplaceOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~LSCutCellLaplaceOperator()

void
LSCutCellLaplaceOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    LS_TIMER_START(t_apply);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<CellVariable<NDIM, double>> x_cc_var = x.getComponentVariable(comp);
        Pointer<CellVariable<NDIM, double>> y_cc_var = y.getComponentVariable(comp);
        if (!x_cc_var || !y_cc_var)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  encountered non-cell centered vector components" << std::endl);
        }
        Pointer<CellDataFactory<NDIM, double>> x_factory = x_cc_var->getPatchDataFactory();
        Pointer<CellDataFactory<NDIM, double>> y_factory = y_cc_var->getPatchDataFactory();
        TBOX_ASSERT(x_factory);
        TBOX_ASSERT(y_factory);
        const unsigned int x_depth = x_factory->getDefaultDepth();
        const unsigned int y_depth = y_factory->getDefaultDepth();
        TBOX_ASSERT(x_depth == y_depth);
        if (x_depth != d_bc_coefs.size() || y_depth != d_bc_coefs.size())
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  each vector component must have data depth == " << d_bc_coefs.size() << "\n"
                                     << "  since d_bc_coefs.size() == " << d_bc_coefs.size() << std::endl);
        }
    }
#endif

    // Loop over comp data
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        {
            using InterpolationTransactionComponent =
                HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<InterpolationTransactionComponent> transaction_comps;
            InterpolationTransactionComponent x_component(d_Q_scr_idx,
                                                          x.getComponentDescriptorIndex(comp),
                                                          DATA_REFINE_TYPE,
                                                          false,
                                                          DATA_COARSEN_TYPE,
                                                          BDRY_EXTRAP_TYPE,
                                                          CONSISTENT_TYPE_2_BDRY,
                                                          d_bc_coefs);
            transaction_comps.push_back(x_component);
            HierarchyGhostCellInterpolation hier_ghost_cell;
            hier_ghost_cell.initializeOperatorState(transaction_comps, d_hierarchy);
            hier_ghost_cell.setHomogeneousBc(d_homogeneous_bc);
            hier_ghost_cell.fillData(d_solution_time);
        }

        Pointer<CellVariable<NDIM, double>> y_cc_var = y.getComponentVariable(comp);
        const int y_idx = y.getComponentDescriptorIndex(comp);
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> x_data = patch->getPatchData(d_Q_scr_idx);
                Pointer<CellData<NDIM, double>> y_data = patch->getPatchData(y_idx);
                computeHelmholtzAction(*x_data, *y_data, *patch);
            }
        }

        if (d_robin_bdry)
        {
            TBOX_ASSERT(d_bdry_conds);
            d_bdry_conds->setHomogeneousBdry(d_homogeneous_bc);
            d_bdry_conds->applyBoundaryCondition(d_Q_var, d_Q_scr_idx, y_cc_var, y_idx, d_hierarchy, d_solution_time);
        }
    }
    LS_TIMER_STOP(t_apply);
    return;
} // apply

void
LSCutCellLaplaceOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                  const SAMRAIVectorReal<NDIM, double>& out)
{
    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in.cloneVector(in.getName());
    d_b = out.cloneVector(out.getName());

    // Setup operator state.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();

    d_ncomp = in.getNumberOfComponents();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
    TBOX_ASSERT(d_ncomp == out.getNumberOfComponents());
#endif

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops =
            new HierarchyMathOps(d_object_name + "::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_hier_math_ops);
#endif
    }

    // Setup the interpolation transaction information.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.clear();
    InterpolationTransactionComponent component(
        d_Q_scr_idx, DATA_REFINE_TYPE, false, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_bc_coefs);
    d_transaction_comps.push_back(component);

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_Q_scr_idx)) level->allocatePatchData(d_Q_scr_idx);
    }

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy, d_coarsest_ln, d_finest_ln);

    if (d_bdry_conds)
    {
        d_bdry_conds->setDiffusionCoefficient(d_poisson_spec.getDConstant());
        d_bdry_conds->setTimeStepType(d_ts_type);
        d_bdry_conds->setLSData(d_ls_var, d_ls_idx, d_vol_var, d_vol_idx, d_area_var, d_area_idx);
        d_bdry_conds->allocateOperatorState(d_hierarchy, d_current_time);
    }

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
} // initializeOperatorState

void
LSCutCellLaplaceOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_Q_scr_idx)) level->deallocatePatchData(d_Q_scr_idx);
    }

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Delete the solution and rhs vectors.
    d_x->freeVectorComponents();
    d_x.setNull();

    d_b->freeVectorComponents();
    d_b.setNull();

    if (d_bdry_conds) d_bdry_conds->deallocateOperatorState(d_hierarchy, d_current_time);

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;
    return;
} // deallocateOperatorState

void
LSCutCellLaplaceOperator::setLSIndices(int ls_idx,
                                       Pointer<NodeVariable<NDIM, double>> ls_var,
                                       int vol_idx,
                                       Pointer<CellVariable<NDIM, double>> vol_var,
                                       int area_idx,
                                       Pointer<CellVariable<NDIM, double>> area_var,
                                       int side_idx,
                                       Pointer<SideVariable<NDIM, double>> side_var)
{
    d_ls_idx = ls_idx;
    d_ls_var = ls_var;
    d_vol_idx = vol_idx;
    d_vol_var = vol_var;
    d_area_idx = area_idx;
    d_area_var = area_var;
    d_side_idx = side_idx;
    d_side_var = side_var;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LSCutCellLaplaceOperator::computeHelmholtzAction(const CellData<NDIM, double>& Q_data,
                                                 CellData<NDIM, double>& R_data,
                                                 const Patch<NDIM>& patch)
{
    LS_TIMER_START(t_compute_helmholtz);
    const Box<NDIM>& box = patch.getBox();
    const Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    const double C = d_poisson_spec.getCConstant();
    const double D =
        d_poisson_spec.dIsConstant() ? d_poisson_spec.getDConstant() : std::numeric_limits<double>::quiet_NaN();

    Pointer<NodeData<NDIM, double>> phi_n_data = patch.getPatchData(d_ls_idx);
    Pointer<SideData<NDIM, double>> D_data =
        d_poisson_spec.dIsConstant() ? nullptr : patch.getPatchData(d_poisson_spec.getDPatchDataId());
    Pointer<CellData<NDIM, double>> area_data = patch.getPatchData(d_area_idx);
    Pointer<CellData<NDIM, double>> vol_data = patch.getPatchData(d_vol_idx);
    Pointer<SideData<NDIM, double>> side_data = patch.getPatchData(d_side_idx);

#if (NDIM == 2)
    IntVector<NDIM> xp(1, 0), yp(0, 1);
#endif
#if (NDIM == 3)
    IntVector<NDIM> xp(1, 0, 0), yp(0, 1, 0), zp(0, 0, 1);
#endif
    for (CellIterator<NDIM> ci(box); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();
        double cell_volume = (*vol_data)(idx);
        for (int d = 0; d < NDIM; ++d) cell_volume *= dx[d];
        if (MathUtilities<double>::equalEps(cell_volume, 0.0))
        {
            for (unsigned int l = 0; l < d_bc_coefs.size(); ++l) R_data(idx, l) = 0.0;
            continue;
        }
        double Q_avg = (Q_data)(idx);
        for (unsigned int l = 0; l < d_bc_coefs.size(); ++l)
        {
            R_data(idx, l) = 0.0;
            // TODO: For NDIM == 3 Update length_fraction to calculate partial face areas instead of
            // lengths.
            // Loop through X faces
            for (int f = 0; f < 2; ++f)
            {
                const int sgn = f == 0 ? -1 : 1;
                const double L = dx[1] * (*side_data)(SideIndex<NDIM>(idx, 0, f));
                double Q_next = Q_data(idx + xp * sgn);
                double vol_next = (*vol_data)(idx + xp * sgn);
                if (MathUtilities<double>::equalEps(vol_next, 0.0)) continue;
                double dudx = (Q_next - Q_avg) / dx[0];
                R_data(idx, l) += D * L * dudx / cell_volume;
            }
            // Loop through Y faces
            for (int f = 0; f < 2; ++f)
            {
                const int sgn = f == 0 ? -1 : 1;
                const double L = dx[0] * (*side_data)(SideIndex<NDIM>(idx, 1, f));
                double Q_next = Q_data(idx + yp * sgn);
                double vol_next = (*vol_data)(idx + yp * sgn);
                if (MathUtilities<double>::equalEps(vol_next, 0.0)) continue;
                double dudy = (Q_next - Q_avg) / dx[1];
                R_data(idx, l) += D * L * dudy / cell_volume;
            }
#if (NDIM == 3)
            for (int f = 0; f < 2; ++f)
            {
                const int sgn = f == 0 ? -1 : 1;
                // Note the awkward ordering of indices. Needs to start at "bottom left" index and go clockwise
                const double L = (*side_data)(SideIndex<NDIM>(idx, 2, f));
                double Q_next = Q_data(idx + zp * sgn);
                double vol_next = (*vol_data)(idx + zp * sgn);
                if (MathUtilities<double>::equalEps(vol_next, 0.0)) continue;
                double dudz = (Q_next - Q_avg) / dx[2];
                R_data(idx, l) += D * L * dudz / cell_volume;
            }
#endif
            // Add C constant
            R_data(idx, l) += C * Q_data(idx, l);
        }
    }
    LS_TIMER_STOP(t_compute_helmholtz);
    return;
}
//////////////////////////////////////////////////////////////////////////////

} // namespace LS

//////////////////////////////////////////////////////////////////////////////
