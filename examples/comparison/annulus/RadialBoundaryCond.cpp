#include "LS/utility_functions.h"

#include "RadialBoundaryCond.h"

namespace
{
static Timer* s_apply_timer = nullptr;
}
RadialBoundaryCond::RadialBoundaryCond(const std::string& object_name, Pointer<Database> input_db)
    : LSCutCellBoundaryConditions(object_name)
{
    input_db->getDoubleArray("Center", d_cent.data(), NDIM);
    d_m = input_db->getDouble("m");
    IBTK_DO_ONCE(s_apply_timer =
                     TimerManager::getManager()->getTimer("LS::RadialBoundaryCond::applyBoundaryCondition"));
}

void
RadialBoundaryCond::applyBoundaryCondition(Pointer<CellVariable<NDIM, double>> Q_var,
                                           const int Q_idx,
                                           Pointer<CellVariable<NDIM, double>> R_var,
                                           const int R_idx,
                                           Pointer<PatchHierarchy<NDIM>> hierarchy,
                                           const double time)
{
    LS_TIMER_START(s_apply_timer);
    TBOX_ASSERT(d_ls_var && d_vol_var && d_area_var);
    TBOX_ASSERT(d_ls_idx > 0 && d_vol_idx > 0 && d_area_idx > 0);

    const double sgn = d_D / std::abs(d_D);
    double pre_fac = sgn * (d_ts_type == LS::DiffusionTimeIntegrationMethod::TRAPEZOIDAL_RULE ? 0.5 : 1.0);
    if (d_D == 0.0) pre_fac = 0.0;

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            const hier::Index<NDIM>& idx_low = box.lower();

            Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(d_area_idx);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
            Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);

            for (CellIterator<NDIM> ci(box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                VectorNd X;
                for (int d = 0; d < NDIM; ++d) X[d] = xlow[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
                const double cell_volume = dx[0] * dx[1] *
#if (NDIM == 3)
                                           dx[2] *
#endif
                                           (*vol_data)(idx);
                const double area = (*area_data)(idx);
                if (area > 0.0)
                {
                    TBOX_ASSERT(cell_volume > 0.0);
                    NodeIndex<NDIM> idx_ll(idx, NodeIndex<NDIM>::LowerLeft);
                    NodeIndex<NDIM> idx_lr(idx, NodeIndex<NDIM>::LowerRight);
                    NodeIndex<NDIM> idx_ul(idx, NodeIndex<NDIM>::UpperLeft);
                    NodeIndex<NDIM> idx_ur(idx, NodeIndex<NDIM>::UpperRight);
                    double dphi_dx =
                        ((*ls_data)(idx_ur) + (*ls_data)(idx_lr) - (*ls_data)(idx_ul) - (*ls_data)(idx_ll)) /
                        (2.0 * dx[0]);
                    double dphi_dy =
                        ((*ls_data)(idx_ul) + (*ls_data)(idx_ur) - (*ls_data)(idx_ll) - (*ls_data)(idx_lr)) /
                        (2.0 * dx[1]);
                    double dist = LS::node_to_cell(idx, *ls_data) / std::sqrt(dphi_dx * dphi_dx + dphi_dy * dphi_dy);
                    double g = d_m * ((X - d_cent).norm() < (0.5 * M_PI) ? -3.0 : 1.0);
                    for (int l = 0; l < Q_data->getDepth(); ++l)
                    {
                        if (!d_homogeneous_bdry)
                        {
                            (*R_data)(idx, l) += pre_fac * g * area / cell_volume;
                        }
                    }
                }
            }
        }
    }
    LS_TIMER_STOP(s_apply_timer);
}
