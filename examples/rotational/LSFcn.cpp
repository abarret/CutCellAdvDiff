// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "LS/utility_functions.h"

#include "LSFcn.h"
/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <SAMRAI_config.h>

#include <array>

namespace LS
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

LSFcn::LSFcn(const string& object_name, Pointer<Database> input_db)
    : CartGridFunction(object_name), d_object_name(object_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif

    // Initialize object with data read from the input database.
    getFromInput(input_db);
    return;
} // LSFcn

LSFcn::~LSFcn()
{
    // intentionally blank
    return;
} // ~LSFcn

void
LSFcn::setDataOnPatch(const int data_idx,
                      Pointer<hier::Variable<NDIM>> /*var*/,
                      Pointer<Patch<NDIM>> patch,
                      const double data_time,
                      const bool /*initial_time*/,
                      Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<CellData<NDIM, double>> ls_c_data = patch->getPatchData(data_idx);
    Pointer<NodeData<NDIM, double>> ls_n_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(ls_c_data || ls_n_data);
#endif
    const Box<NDIM>& patch_box = d_extended_box ? ls_n_data->getGhostBox() : patch->getBox();
    const hier::Index<NDIM>& patch_lower = patch->getBox().lower();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();

    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    VectorNd X;
    // Rotate center of mass to current configuration
    MatrixNd Q;
    Q(0, 0) = Q(1, 1) = std::cos(data_time * (2.0 * M_PI));
    Q(0, 1) = -std::sin(data_time * (2.0 * M_PI));
    Q(1, 0) = std::sin(data_time * (2.0 * M_PI));
    VectorNd current_com = Q * d_com;

    if (ls_c_data)
    {
        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            for (int d = 0; d < NDIM; ++d)
                X[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - patch_lower(d)) + 0.5);
            const double r = (X - current_com).norm();
            (*ls_c_data)(idx) = r - d_R;
        }
    }
    else if (ls_n_data)
    {
        for (NodeIterator<NDIM> ni(patch_box); ni; ni++)
        {
            const NodeIndex<NDIM>& idx = ni();
            for (int d = 0; d < NDIM; ++d) X[d] = x_lower[d] + dx[d] * static_cast<double>(idx(d) - patch_lower(d));
            const double r = (X - current_com).norm();
            (*ls_n_data)(idx) = r - d_R;
        }
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LSFcn::getFromInput(Pointer<Database> db)
{
    db->getDoubleArray("com", d_com.data(), NDIM);
    d_R = db->getDouble("r");
    return;
} // getFromInput
} // namespace LS
