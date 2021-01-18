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

#include "OutsideLSFcn.h"
/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <SAMRAI_config.h>

#include <array>

/////////////////////////////// PUBLIC ///////////////////////////////////////

OutsideLSFcn::OutsideLSFcn(const string& object_name,
                           Pointer<HierarchyIntegrator> hierarchy_integrator,
                           Pointer<NodeVariable<NDIM, double>> in_var,
                           Pointer<Database> input_db)
    : CartGridFunction(object_name), d_integrator(hierarchy_integrator), d_in_ls_node_var(in_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif

    // Initialize object with data read from the input database.
    d_R1 = input_db->getDouble("R1");
    d_R2 = input_db->getDouble("R2");
    return;
} // OutsideLSFcn

OutsideLSFcn::OutsideLSFcn(Pointer<VariableContext> ctx,
                           const std::string& object_name,
                           Pointer<NodeVariable<NDIM, double>> in_n_var,
                           Pointer<Database> input_db)
    : CartGridFunction(object_name),
      d_ctx(ctx),
      d_in_ls_node_var(in_n_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif

    // Initialize object with data read from the input database.
    d_R1 = input_db->getDouble("R1");
    d_R2 = input_db->getDouble("R2");
    return;
}

OutsideLSFcn::~OutsideLSFcn()
{
    // intentionally destroy hierarchy_integrator to prevent circular references
    d_integrator.setNull();
    d_ctx.setNull();
    return;
} // ~OutsideLSFcn

void
OutsideLSFcn::setDataOnPatch(const int data_idx,
                             Pointer<hier::Variable<NDIM>> /*var*/,
                             Pointer<Patch<NDIM>> patch,
                             const double data_time,
                             const bool /*initial_time*/,
                             Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<NodeData<NDIM, double>> ls_out_node_data = patch->getPatchData(data_idx);

    Pointer<VariableContext> ctx;
    if (d_integrator)
    {
        ctx = MathUtilities<double>::equalEps(data_time, d_integrator->getIntegratorTime()) ?
                  d_integrator->getCurrentContext() :
                  d_integrator->getNewContext();
    }
    else if (d_ctx)
    {
        ctx = d_ctx;
    }
    else
    {
        TBOX_ERROR("SHOULD NOT REACH HERE.\n");
    }

    Pointer<NodeData<NDIM, double>> ls_in_node_data =
        d_in_ls_node_var ? patch->getPatchData(d_in_ls_node_var, ctx) : nullptr;

    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_low = pgeom->getXLower();

    const Box<NDIM>& box = patch->getBox();
    const hier::Index<NDIM>& idx_low = box.lower();

    if (ls_out_node_data && ls_in_node_data)
    {
        for (NodeIterator<NDIM> ci(box); ci; ci++)
        {
            const NodeIndex<NDIM>& idx = ci();

            const double ls_in = (*ls_in_node_data)(idx);

            VectorNd x_pt;
            for (int d = 0; d < NDIM; ++d) x_pt(d) = x_low[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
            const double r = x_pt.norm();
            const double ls_disk = std::max(r - d_R2, d_R1 - r);
            (*ls_out_node_data)(idx) = std::max(ls_disk, -ls_in);
        }
    }
    else
    {
        TBOX_ERROR("Should not get here.");
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
