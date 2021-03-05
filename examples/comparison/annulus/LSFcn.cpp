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
    d_R1 = input_db->getDouble("R1");
    d_R2 = input_db->getDouble("R2");
    input_db->getDoubleArray("center", d_cent.data(), NDIM);
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
    Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(data_idx);
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_low = pgeom->getXLower();

    const Box<NDIM>& box = patch->getBox();
    const hier::Index<NDIM>& idx_low = box.lower();
    for (NodeIterator<NDIM> ci(box); ci; ci++)
    {
        const NodeIndex<NDIM>& idx = ci();

        VectorNd x_pt;
        for (int d = 0; d < NDIM; ++d) x_pt(d) = x_low[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
        const double r = (x_pt - d_cent).norm();
        (*ls_data)(idx) = std::max(r - d_R2, d_R1 - r);
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
