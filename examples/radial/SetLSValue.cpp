#include "ibamr/config.h"

#include <ibamr/app_namespaces.h>

#include "SetLSValue.h"

#include <SAMRAI_config.h>

#include <array>

/////////////////////////////// PUBLIC ///////////////////////////////////////

SetLSValue::SetLSValue(const string& object_name, Pointer<GridGeometry<NDIM>> grid_geom, Pointer<Database> input_db)
    : CartGridFunction(object_name), d_grid_geom(grid_geom)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
    TBOX_ASSERT(d_grid_geom);
#endif

    // Initialize object with data read from the input database.
    getFromInput(input_db);
    return;
} // SetLSValue

SetLSValue::~SetLSValue()
{
    // intentionally blank
    return;
} // ~SetLSValue

void
SetLSValue::setDataOnPatch(const int data_idx,
                           Pointer<hier::Variable<NDIM>> /*var*/,
                           Pointer<Patch<NDIM>> patch,
                           const double data_time,
                           const bool /*initial_time*/,
                           Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<NodeData<NDIM, double>> ls_n_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(ls_n_data);
#endif
    const Box<NDIM>& patch_box = d_extended_box ? ls_n_data->getGhostBox() : patch->getBox();
    const hier::Index<NDIM>& patch_lower = patch->getBox().lower();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();

    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    VectorNd X;
    for (NodeIterator<NDIM> ni(patch_box); ni; ni++)
    {
        const NodeIndex<NDIM>& i = ni();
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = x_lower[d] + dx[d] * static_cast<double>(i(d) - patch_lower(d)) - data_time * d_U[d];
        }
        X = X - d_center;
        double r = X.norm();
        (*ls_n_data)(i) = r - d_R1;
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SetLSValue::getFromInput(Pointer<Database> db)
{
    if (db)
    {
        d_R1 = db->getDoubleWithDefault("R1", d_R1);
        if (db->keyExists("Center")) db->getDoubleArray("Center", d_center.data(), NDIM);
        if (db->keyExists("Vel")) db->getDoubleArray("Vel", d_U.data(), NDIM);
    }
    return;
} // getFromInput
