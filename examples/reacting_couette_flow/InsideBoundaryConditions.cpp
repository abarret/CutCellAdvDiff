#include "InsideBoundaryConditions.h"

#include "ibamr/app_namespaces.h"

InsideBoundaryConditions::InsideBoundaryConditions(const std::string& object_name,
                                                   Pointer<Database> input_db,
                                                   Pointer<CellVariable<NDIM, double>> out_var,
                                                   Pointer<AdvDiffHierarchyIntegrator> integrator)
    : LSCutCellBoundaryConditions(object_name), d_out_var(out_var), d_integrator(integrator)
{
    d_k1 = input_db->getDouble("k1");
}

InsideBoundaryConditions::~InsideBoundaryConditions()
{
    // automatically deallocate cell variable and integrator to prevent circular definitions
    d_out_var.setNull();
    d_integrator.setNull();
}

void
InsideBoundaryConditions::applyBoundaryCondition(Pointer<CellVariable<NDIM, double>> Q_var,
                                                 const int Q_idx,
                                                 Pointer<CellVariable<NDIM, double>> R_var,
                                                 const int R_idx,
                                                 Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                 const double time)
{
    TBOX_ASSERT(d_ls_var && d_vol_var && d_area_var);
    TBOX_ASSERT(d_ls_idx > 0 && d_vol_idx > 0 && d_area_idx > 0);

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int out_idx = var_db->mapVariableAndContextToIndex(d_out_var, d_integrator->getCurrentContext());

    const double sgn = d_D / std::abs(d_D);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            Pointer<CellData<NDIM, double>> in_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double>> out_data = patch->getPatchData(out_idx);
            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(d_area_idx);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);

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
                    if (!d_homogeneous_bdry) (*R_data)(idx) += 0.5 * sgn * d_k1 * (*out_data)(idx)*area / cell_volume;
                    (*R_data)(idx) -= 0.5 * sgn * d_k1 * (*in_data)(idx)*area / cell_volume;
                }
            }
        }
    }
}
