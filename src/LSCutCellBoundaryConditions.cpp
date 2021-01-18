#include "LS/LSCutCellBoundaryConditions.h"

#include "ibamr/app_namespaces.h"

namespace LS
{
LSCutCellBoundaryConditions::LSCutCellBoundaryConditions(const std::string& object_name) : d_object_name(object_name)
{
    // intentionally blank
}

void
LSCutCellBoundaryConditions::allocateOperatorState(Pointer<PatchHierarchy<NDIM>> /*hierarchy*/, double /*time*/)
{
    // intentionally blank
}

void
LSCutCellBoundaryConditions::deallocateOperatorState(Pointer<PatchHierarchy<NDIM>> /*hierarchy*/, double /*time*/)
{
    // intentionally blank
}

void
LSCutCellBoundaryConditions::setLSData(Pointer<NodeVariable<NDIM, double>> ls_var,
                                       const int ls_idx,
                                       Pointer<CellVariable<NDIM, double>> vol_var,
                                       const int vol_idx,
                                       Pointer<CellVariable<NDIM, double>> area_var,
                                       const int area_idx)
{
    d_ls_var = ls_var;
    d_ls_idx = ls_idx;
    d_vol_var = vol_var;
    d_vol_idx = vol_idx;
    d_area_var = area_var;
    d_area_idx = area_idx;
}

void
LSCutCellBoundaryConditions::setHomogeneousBdry(const bool homogeneous_bdry)
{
    d_homogeneous_bdry = homogeneous_bdry;
}

void
LSCutCellBoundaryConditions::setDiffusionCoefficient(const double D)
{
    d_D = D;
}
} // namespace LS
