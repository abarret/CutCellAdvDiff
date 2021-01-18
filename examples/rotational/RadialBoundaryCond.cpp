#include "RadialBoundaryCond.h"

#include "ibamr/app_namespaces.h"

RadialBoundaryCond::RadialBoundaryCond(const std::string& object_name, Pointer<Database> input_db)
    : LSCutCellBoundaryConditions(object_name)
{
    d_a = input_db->getDouble("a_coef");
    d_D_coef = input_db->getDouble("D");
    d_R = input_db->getDouble("R");
}

void
RadialBoundaryCond::applyBoundaryCondition(Pointer<CellVariable<NDIM, double>> Q_var,
                                           const int Q_idx,
                                           Pointer<CellVariable<NDIM, double>> R_var,
                                           const int R_idx,
                                           Pointer<PatchHierarchy<NDIM>> hierarchy,
                                           const double time)
{
    TBOX_ASSERT(d_ls_var && d_vol_var && d_area_var);
    TBOX_ASSERT(d_ls_idx > 0 && d_vol_idx > 0 && d_area_idx > 0);

    double g = 5.0 * (d_a - d_R + 2.0 * d_a * time) /
               (d_D_coef * std::exp(d_R * d_R / (2.0 * d_D_coef + 4.0 * d_D_coef * time)) * (1.0 + 2.0 * time) *
                (1.0 + 2.0 * time));

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
                    for (int l = 0; l < Q_data->getDepth(); ++l)
                    {
                        if (d_homogeneous_bdry)
                        {
                            (*R_data)(idx, l) -= pre_fac * d_a * area *
                                                 ((*Q_data)(idx, l) * d_D_coef / (d_D_coef - d_a * dist)) / cell_volume;
                        }
                        else
                        {
                            (*R_data)(idx, l) += 0.5 * sgn * g * area / cell_volume;
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
}
