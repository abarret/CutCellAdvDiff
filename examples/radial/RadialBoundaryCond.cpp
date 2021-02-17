#include "LS/utility_functions.h"

#include "RadialBoundaryCond.h"

namespace
{
static Timer* s_apply_timer = nullptr;
}
RadialBoundaryCond::RadialBoundaryCond(const std::string& object_name, Pointer<Database> input_db)
    : LSCutCellBoundaryConditions(object_name)
{
    d_a = input_db->getDouble("a_coef");
    d_D_coef = input_db->getDouble("D");
    d_R = input_db->getDouble("R");
    input_db->getDoubleArray("Center", d_center.data(), NDIM);
    d_vel.resize(NDIM);
    input_db->getDoubleArray("velocity", d_vel.data(), NDIM);
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
#if (NDIM == 2)
    double g = 5.0 * (d_a - d_R + 2.0 * d_a * time) /
               (d_D_coef * std::exp(d_R * d_R / (2.0 * d_D_coef + 4.0 * d_D_coef * time)) * M_PI * (1.0 + 2.0 * time) *
                (1.0 + 2.0 * time));
#endif
#if (NDIM == 3)
    double g = 5.0 * (d_a - d_R + 2.0 * d_a * time) /
               (d_D_coef * std::exp(d_R * d_R / (2.0 * d_D_coef + 4.0 * d_D_coef * time)) * (1.0 + 2.0 * time) *
                (1.0 + 2.0 * time) * std::sqrt(2.0 * M_PI * M_PI * M_PI * (d_D_coef + 2.0 * d_D_coef * time)));
#endif

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

            Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(d_area_idx);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
            Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);

            for (CellIterator<NDIM> ci(box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                const double cell_volume = dx[0] * dx[1] *
#if (NDIM == 3)
                                           dx[2] *
#endif
                                           (*vol_data)(idx);
                const double area = (*area_data)(idx);
                if (area > 0.0)
                {
                    TBOX_ASSERT(cell_volume > 0.0);
#if (NDIM == 2)
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
#endif
#if (NDIM == 3)
                    NodeIndex<NDIM> idx_lll(idx, NodeIndex<NDIM>::LLL);
                    NodeIndex<NDIM> idx_llu(idx, NodeIndex<NDIM>::LLU);
                    NodeIndex<NDIM> idx_lul(idx, NodeIndex<NDIM>::LUL);
                    NodeIndex<NDIM> idx_luu(idx, NodeIndex<NDIM>::LUU);
                    NodeIndex<NDIM> idx_ull(idx, NodeIndex<NDIM>::ULL);
                    NodeIndex<NDIM> idx_ulu(idx, NodeIndex<NDIM>::ULU);
                    NodeIndex<NDIM> idx_uul(idx, NodeIndex<NDIM>::UUL);
                    NodeIndex<NDIM> idx_uuu(idx, NodeIndex<NDIM>::UUU);
                    VectorNd grad_phi;
                    grad_phi[0] =
                        (((*ls_data)(idx_ull)) + (*ls_data)(idx_ulu) + (*ls_data)(idx_uul) + (*ls_data)(idx_uuu) -
                         ((*ls_data)(idx_lll) + (*ls_data)(idx_llu) + (*ls_data)(idx_lul) + (*ls_data)(idx_luu))) /
                        (0.25 * dx[0]);
                    grad_phi[1] =
                        (((*ls_data)(idx_lul)) + (*ls_data)(idx_luu) + (*ls_data)(idx_uul) + (*ls_data)(idx_uuu) -
                         ((*ls_data)(idx_lll) + (*ls_data)(idx_llu) + (*ls_data)(idx_ull) + (*ls_data)(idx_ulu))) /
                        (0.25 * dx[1]);
                    grad_phi[2] =
                        (((*ls_data)(idx_llu)) + (*ls_data)(idx_luu) + (*ls_data)(idx_ulu) + (*ls_data)(idx_uuu) -
                         ((*ls_data)(idx_lll) + (*ls_data)(idx_lul) + (*ls_data)(idx_ull) + (*ls_data)(idx_uul))) /
                        (0.25 * dx[2]);
                    double dist = LS::node_to_cell(idx, *ls_data) / grad_phi.norm();
#endif
                    for (int l = 0; l < Q_data->getDepth(); ++l)
                    {
                        if (d_homogeneous_bdry)
                        {
                            (*R_data)(idx, l) -= pre_fac * d_a * area *
                                                 ((*Q_data)(idx, l) * d_D_coef / (d_D_coef - d_a * dist)) / cell_volume;
                        }
                        else
                        {
                            (*R_data)(idx, l) += pre_fac * g * area / cell_volume;
                            (*R_data)(idx, l) -= pre_fac * d_a * area *
                                                 ((*Q_data)(idx, l) * d_D_coef / (d_D_coef - d_a * dist) -
                                                  g * dist / (d_D_coef - d_a * dist)) /
                                                 cell_volume;
                        }
                    }
                }
            }
        }
    }
    LS_TIMER_STOP(s_apply_timer);
}
