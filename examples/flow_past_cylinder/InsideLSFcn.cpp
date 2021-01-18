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

#include <ibamr/config.h>

#include "LS/utility_functions.h"

#include "InsideLSFcn.h"

#include <SAMRAI_config.h>

#include <array>

#include "ibamr/app_namespaces.h"

double InsideLSFcn::s_large_num = 0.5;

namespace
{
static Timer* t_fill_ls;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

InsideLSFcn::InsideLSFcn(const string& object_name,
                         Pointer<Database> input_db,
                         BoundaryMesh* ls_mesh,
                         BoundaryMesh* disk_mesh,
                         FEDataManager* disk_fe_manager)
    : CartGridFunction(object_name),
      d_ls_scr_var(new CellVariable<NDIM, double>(d_object_name + "::LS_SCR")),
      d_ls_mesh(ls_mesh),
      d_disk_mesh(disk_mesh),
      d_disk_fe_manager(disk_fe_manager)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif
    input_db->getDoubleArray("center", d_cent.data(), NDIM);
    d_R = input_db->getDouble("radius");

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_ls_scr_idx =
        var_db->registerVariableAndContext(d_ls_scr_var, var_db->getContext(d_object_name + "::Distance"), 1);
    d_n_scr_idx = var_db->registerVariableAndContext(d_ls_scr_var, var_db->getContext(d_object_name + "::num_elems"));
    IBAMR_DO_ONCE(t_fill_ls = TimerManager::getManager()->getTimer("LS::InsideLSFcn::setDataOnPatchHierarchy"););
    return;
} // InsideLSFcn

InsideLSFcn::~InsideLSFcn()
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_ls_scr_idx)) level->deallocatePatchData(d_ls_scr_idx);
    }
}

void
InsideLSFcn::setDataOnPatchHierarchy(const int data_idx,
                                     Pointer<hier::Variable<NDIM>> var,
                                     Pointer<PatchHierarchy<NDIM>> hierarchy,
                                     const double data_time,
                                     const bool initial_time,
                                     int coarsest_ln,
                                     int finest_ln)
{
    LS_TIMER_START(t_fill_ls);
    d_hierarchy = hierarchy;
    coarsest_ln = coarsest_ln == -1 ? 0 : coarsest_ln;
    finest_ln = finest_ln == -1 ? hierarchy->getFinestLevelNumber() : finest_ln;

    // If we cached time is not data_time, we need to update the mesh, and recreate the level set function
    if (d_cached_time != data_time)
    {
        d_cached_time = data_time;
        if (!d_surface_distance_evaluator)
            d_surface_distance_evaluator = new FESurfaceDistanceEvaluator(
                d_object_name + "::FESurfaceDistanceEvaluator", hierarchy, *d_ls_mesh, *d_ls_mesh, 2, true);
        updateMesh(*d_ls_mesh, *d_disk_mesh, d_disk_fe_manager);
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(d_ls_scr_idx)) level->allocatePatchData(d_ls_scr_idx);
            if (!level->checkAllocated(d_n_scr_idx)) level->allocatePatchData(d_n_scr_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> ls_data = patch->getPatchData(d_ls_scr_idx);
                ls_data->fillAll(s_large_num);
            }
        }
        pout << "Started mapping intersections" << std::endl;
        d_surface_distance_evaluator->mapIntersections();
        pout << "Finished mapping intersections" << std::endl;
        pout << "Computing face normal" << std::endl;
        d_surface_distance_evaluator->calculateSurfaceNormals();
        pout << "Finished calculation of face normal" << std::endl;
        pout << "Computing distances" << std::endl;
        d_surface_distance_evaluator->computeSignedDistance(d_n_scr_idx, d_ls_scr_idx);
        pout << "Finished computing distances" << std::endl;
        pout << "Updating sign away from interface" << std::endl;
        d_surface_distance_evaluator->updateSignAwayFromInterface(d_ls_scr_idx, d_hierarchy, s_large_num);
        pout << "Finished updating sign" << std::endl;
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comps(1);
        ghost_cell_comps[0] = ITC(d_ls_scr_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR");
        HierarchyGhostCellInterpolation hier_ghost_cells;
        hier_ghost_cells.initializeOperatorState(ghost_cell_comps, d_hierarchy);
        hier_ghost_cells.fillData(data_time);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        setDataOnPatchLevel(data_idx, var, level, data_time, initial_time);
    }

    LS_TIMER_STOP(t_fill_ls);
}

void
InsideLSFcn::setDataOnPatch(const int data_idx,
                            Pointer<hier::Variable<NDIM>> var,
                            Pointer<Patch<NDIM>> patch,
                            const double data_time,
                            const bool initial_time,
                            Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<PatchData<NDIM>> data = patch->getPatchData(data_idx);

    if (initial_time)
    {
        initialLSValue(data, var, patch, data_time);
    }
    else
    {
        Pointer<CellData<NDIM, double>> ls_data = patch->getPatchData(d_ls_scr_idx);
        Pointer<CellData<NDIM, double>> c_data = data;
        Pointer<NodeData<NDIM, double>> n_data = data;
        if (n_data)
        {
            // Interpolate cell data to nodes
            PatchMathOps patch_ops;
            if (ls_data) patch_ops.interp(n_data, ls_data, patch, false);
            // Now negate the data, since we want the outside disk
            for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
            {
                const NodeIndex<NDIM>& idx = ni();
                (*n_data)(idx) = -(*n_data)(idx);
            }
        }
        else if (c_data)
        {
            // Copy cell data
            if (ls_data) c_data->copy(*ls_data);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                (*c_data)(idx) = -(*c_data)(idx);
            }
        }
        else
        {
            TBOX_ERROR("Should not get here.\n");
        }
    }
} // setDataOnPatch

void
InsideLSFcn::updateMesh(Mesh& ls_mesh, Mesh& disk_mesh, FEDataManager* disk_fe_manager)
{
    plog << "Updating mesh.\n";
    EquationSystems* eq_sys = disk_fe_manager->getEquationSystems();
    System& X_sys = eq_sys->get_system(disk_fe_manager->COORDINATES_SYSTEM_NAME);
    DofMap& X_dof_map = X_sys.get_dof_map();
    NumericVector<double>* X_vec = X_sys.solution.get();
    // Loop through lower and upper mesh
    const MeshBase::const_node_iterator end_n = ls_mesh.nodes_end();
    for (MeshBase::const_node_iterator nl = ls_mesh.nodes_begin(); nl != end_n; ++nl)
    {
        Node* const n = *nl;
        Node* const disk_n = disk_mesh.node_ptr(n->id());
        std::vector<dof_id_type> X_dofs;
        for (int d = 0; d < NDIM; ++d)
        {
            IBTK::get_nodal_dof_indices(X_dof_map, disk_n, d, X_dofs);
            double X_val;
            X_vec->get(X_dofs, &X_val);
            (*n)(d) = X_val;
        }
    }
    X_vec->close();
    ls_mesh.prepare_for_use();
    plog << "Finished updating mesh.\n";
}

void
InsideLSFcn::initialLSValue(Pointer<PatchData<NDIM>> data,
                            Pointer<hier::Variable<NDIM>> /*var*/,
                            Pointer<Patch<NDIM>> patch,
                            const double /*data_time*/)
{
    Pointer<CellData<NDIM, double>> c_data = data;
    Pointer<NodeData<NDIM, double>> n_data = data;
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_low = pgeom->getXLower();
    const Box<NDIM>& box = patch->getBox();
    const hier::Index<NDIM>& idx_low = box.lower();

    if (c_data)
    {
        for (CellIterator<NDIM> ci(box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();

            VectorNd x_pt;
            for (int d = 0; d < NDIM; ++d)
                x_pt(d) = x_low[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
            (*c_data)(idx) = d_R - (x_pt - d_cent).norm();
        }
    }
    else if (n_data)
    {
        for (NodeIterator<NDIM> ci(box); ci; ci++)
        {
            const NodeIndex<NDIM>& idx = ci();

            VectorNd x_pt;
            for (int d = 0; d < NDIM; ++d) x_pt(d) = x_low[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
            (*n_data)(idx) = d_R - (x_pt - d_cent).norm();
        }
    }
    else
    {
        TBOX_ERROR("Should not get here.");
    }
}
//////////////////////////////////////////////////////////////////////////////
