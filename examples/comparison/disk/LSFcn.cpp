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

#include <IBAMR_config.h>

#include <SAMRAI_config.h>

#include <array>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LSFcn::LSFcn(const string& object_name, Pointer<Database> input_db) : CartGridFunction(object_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif

    // Initialize object with data read from the input database.
    d_R = input_db->getDouble("r");
    input_db->getDoubleArray("x_cent", d_x_cent.data(), NDIM);
    return;
} // LSFcn

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

    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_low = pgeom->getXLower();

    const Box<NDIM>& box = patch->getBox();
    const hier::Index<NDIM>& idx_low = box.lower();

    if (ls_c_data)
    {
        for (CellIterator<NDIM> ci(box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();

            VectorNd x_pt;
            for (int d = 0; d < NDIM; ++d)
                x_pt(d) = x_low[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
            x_pt -= d_x_cent;
            (*ls_c_data)(idx) = d_R - (x_pt - d_x_cent).norm();
        }
    }
    else if (ls_n_data)
    {
        for (NodeIterator<NDIM> ci(box); ci; ci++)
        {
            const NodeIndex<NDIM>& idx = ci();

            VectorNd x_pt;
            for (int d = 0; d < NDIM; ++d) x_pt(d) = x_low[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
            x_pt -= d_x_cent;
            (*ls_n_data)(idx) = d_R - (x_pt - d_x_cent).norm();
        }
    }
    else
    {
        TBOX_ERROR("Should not get here.");
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
