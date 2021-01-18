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

#include "InsideLSFcn.h"
/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <SAMRAI_config.h>

#include <array>

/////////////////////////////// PUBLIC ///////////////////////////////////////

InsideLSFcn::InsideLSFcn(const string& object_name, Pointer<Database> input_db)
    : CartGridFunction(object_name), d_object_name(object_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif

    // Initialize object with data read from the input database.
    d_R = input_db->getDouble("r");
    d_a = input_db->getDouble("a");
    d_b = input_db->getDouble("b");
    d_period = input_db->getDouble("period");
    input_db->getDoubleArray("x_cent", d_x_cent.data(), NDIM);
    return;
} // InsideLSFcn

void
InsideLSFcn::setDataOnPatch(const int data_idx,
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
            const double r = x_pt.norm() + 1.0e-12;
            const double th = std::atan2(static_cast<double>(x_pt[1]), static_cast<double>(x_pt[0]));
            double new_th = (d_a * r + d_b / r) * (d_period / (2.0 * M_PI)) * cos(2.0 * M_PI * data_time / d_period) +
                            (th - (d_a * r + d_b / r) * (d_period / (2.0 * M_PI)));
            x_pt[0] = r * std::cos(new_th);
            x_pt[1] = r * std::sin(new_th);
            (*ls_c_data)(idx) = (x_pt - d_x_cent).norm() - d_R;
        }
    }
    else if (ls_n_data)
    {
        for (NodeIterator<NDIM> ci(box); ci; ci++)
        {
            const NodeIndex<NDIM>& idx = ci();

            VectorNd x_pt;
            for (int d = 0; d < NDIM; ++d) x_pt(d) = x_low[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
            const double r = x_pt.norm() + 1.0e-12;
            const double th = std::atan2(static_cast<double>(x_pt[1]), static_cast<double>(x_pt[0]));
            double new_th = (d_a * r + d_b / r) * (d_period / (2.0 * M_PI)) * cos(2.0 * M_PI * data_time / d_period) +
                            (th - (d_a * r + d_b / r) * (d_period / (2.0 * M_PI)));
            x_pt[0] = r * std::cos(new_th);
            x_pt[1] = r * std::sin(new_th);
            (*ls_n_data)(idx) = (x_pt - d_x_cent).norm() - d_R;
        }
    }
    else
    {
        TBOX_ERROR("Should not get here.");
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
