#include "ibtk/IndexUtilities.h"

#include "LS/SBBoundaryConditions.h"
#include "LS/utility_functions.h"

#include "libmesh/elem_cutter.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/transient_system.h"

#include <boost/multi_array.hpp>

namespace
{
static Timer* t_applyBoundaryCondition = nullptr;
static Timer* t_interpolateToBoundary = nullptr;
static Timer* t_allocateOperatorState = nullptr;
static Timer* t_deallocateOperatorState = nullptr;
static Timer* t_cache_boundary_data = nullptr;
} // namespace

namespace LS
{
SBBoundaryConditions::SBBoundaryConditions(const std::string& object_name,
                                           Pointer<Database> input_db,
                                           Mesh* mesh,
                                           FEDataManager* fe_data_manager,
                                           bool use_old_soln)
    : LSCutCellBoundaryConditions(object_name),
      d_mesh(mesh),
      d_fe_data_manager(fe_data_manager),
      d_sys_name(d_object_name),
      d_scr_var(new CellVariable<NDIM, double>(d_object_name + "::SCR")),
      d_use_old_soln(use_old_soln)
{
    d_perturb_nodes = input_db->getBool("perturb_nodes");
    ExplicitSystem& fl_sys = d_fe_data_manager->getEquationSystems()->add_system<ExplicitSystem>(d_sys_name);
    fl_sys.add_variable(d_object_name + "::scratch");
    fl_sys.assemble_before_solve = false;
    fl_sys.assemble();

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_scr_idx =
        var_db->registerVariableAndContext(d_scr_var, var_db->getContext(d_object_name + "::SCR"), IntVector<NDIM>(2));

    IBTK_DO_ONCE(t_applyBoundaryCondition =
                     TimerManager::getManager()->getTimer("LS::SBBoundaryConditions::applyBoundaryCondition()");
                 t_interpolateToBoundary =
                     TimerManager::getManager()->getTimer("LS::SBBoundaryConditions::interpolateToBoundary()");
                 t_allocateOperatorState =
                     TimerManager::getManager()->getTimer("LS::SBBoundaryConditions::allocateOperatorState()");
                 t_deallocateOperatorState =
                     TimerManager::getManager()->getTimer("LS::SBBoundaryConditions::deallocateOperatorState()");
                 t_cache_boundary_data =
                     TimerManager::getManager()->getTimer("LS::SBBoundaryConditions::cacheBoundaryData()"););
}

void
SBBoundaryConditions::registerFluidFluidInteraction(Pointer<CellVariable<NDIM, double>> fl_var)
{
    TBOX_ASSERT(fl_var);
    TBOX_ASSERT(std::find(d_fl_vars.begin(), d_fl_vars.end(), fl_var) == d_fl_vars.end());
    d_fl_vars.push_back(fl_var);
    std::string fl_sys_name = d_object_name + "::" + fl_var->getName();
    d_fl_names.push_back(fl_sys_name);
    ExplicitSystem& fl_sys = d_fe_data_manager->getEquationSystems()->add_system<ExplicitSystem>(fl_sys_name);
    fl_sys.add_variable(fl_var->getName());
    fl_sys.assemble_before_solve = false;
    fl_sys.assemble();
}

void
SBBoundaryConditions::setFluidContext(Pointer<VariableContext> ctx)
{
    d_ctx = ctx;
}

void
SBBoundaryConditions::registerFluidSurfaceInteraction(const std::string& surface_name)
{
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), surface_name) == d_sf_names.end());
    d_sf_names.push_back(surface_name);
}

void
SBBoundaryConditions::allocateOperatorState(Pointer<PatchHierarchy<NDIM>> hierarchy, double time)
{
    LS_TIMER_START(t_allocateOperatorState);
    LSCutCellBoundaryConditions::allocateOperatorState(hierarchy, time);
    d_hierarchy = hierarchy;

    TBOX_ASSERT(d_ctx);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_scr_idx)) level->allocatePatchData(d_scr_idx);
    }
    // Interpolate to boundary
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_rbf_reconstruct.setLSData(d_ls_idx, d_vol_idx);
    d_rbf_reconstruct.setPatchHierarchy(hierarchy);
    d_rbf_reconstruct.setUseCentroids(false);

    cacheBoundaryData();

    for (size_t l = 0; l < d_fl_names.size(); ++l)
    {
        const int fl_idx = var_db->mapVariableAndContextToIndex(d_fl_vars[l], d_ctx);
        interpolateToBoundary(fl_idx, d_fl_names[l], hierarchy, time);
    }
    LS_TIMER_STOP(t_allocateOperatorState);
}

void
SBBoundaryConditions::deallocateOperatorState(Pointer<PatchHierarchy<NDIM>> hierarchy, double time)
{
    LS_TIMER_START(t_deallocateOperatorState);
    LSCutCellBoundaryConditions::deallocateOperatorState(hierarchy, time);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_scr_idx)) level->deallocatePatchData(d_scr_idx);
    }
    d_rbf_reconstruct.clearCache();
    clearCache();
    LS_TIMER_STOP(t_deallocateOperatorState);
}

void
SBBoundaryConditions::applyBoundaryCondition(Pointer<CellVariable<NDIM, double>> Q_var,
                                             const int Q_idx,
                                             Pointer<CellVariable<NDIM, double>> R_var,
                                             const int R_idx,
                                             Pointer<PatchHierarchy<NDIM>> hierarchy,
                                             const double time)
{
    if (!d_is_allocated) TBOX_ERROR(d_object_name + "::applyBoundaryCondition() must allocate operator first.\n\n");
    LS_TIMER_START(t_applyBoundaryCondition);
    TBOX_ASSERT(d_ls_var && d_vol_var && d_area_var);
    TBOX_ASSERT(d_ls_idx > 0 && d_vol_idx > 0 && d_area_idx > 0);

    interpolateToBoundary(Q_idx, d_sys_name, hierarchy, time);

    const double sgn = d_D / std::abs(d_D);
    double pre_fac = sgn * (d_ts_type == LS::DiffusionTimeIntegrationMethod::TRAPEZOIDAL_RULE ? 0.5 : 1.0);
    if (d_D == 0.0) pre_fac = 0.0;

    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();

    System& X_system = eq_sys->get_system(d_fe_data_manager->COORDINATES_SYSTEM_NAME);
    DofMap& X_dof_map = X_system.get_dof_map();
    FEType X_fe_type = X_dof_map.variable_type(0);
    NumericVector<double>* X_vec = X_system.solution.get();
    auto X_petsc_vec = dynamic_cast<PetscVector<double>*>(X_vec);
    TBOX_ASSERT(X_petsc_vec != nullptr);
    const double* const X_local_soln = X_petsc_vec->get_array_read();
    FEDataManager::SystemDofMapCache& X_dof_map_cache =
        *d_fe_data_manager->getDofMapCache(d_fe_data_manager->COORDINATES_SYSTEM_NAME);

    std::vector<NumericVector<double>*> fl_vecs, sf_vecs;
    std::vector<DofMap*> fl_dof_maps, sf_dof_maps;
    for (const auto& fl_name : d_fl_names)
    {
        System& fl_sys = eq_sys->get_system(fl_name);
        fl_dof_maps.push_back(&fl_sys.get_dof_map());
        fl_vecs.push_back(fl_sys.solution.get());
    }

    for (const auto sf_name : d_sf_names)
    {
        auto& sf_sys = eq_sys->get_system<TransientExplicitSystem>(sf_name);
        sf_dof_maps.push_back(&sf_sys.get_dof_map());
        if (d_use_old_soln)
            sf_vecs.push_back(sf_sys.old_local_solution.get());
        else
            sf_vecs.push_back(sf_sys.current_local_solution.get());
    }

    // Get base system
    System& Q_sys = eq_sys->get_system(d_sys_name);
    DofMap& Q_dof_map = Q_sys.get_dof_map();
    FEType Q_fe_type = Q_dof_map.variable_type(0);
    NumericVector<double>* Q_vec = Q_sys.solution.get();
    TBOX_ASSERT(Q_fe_type == X_fe_type);

    std::unique_ptr<FEBase> fe = FEBase::build(d_mesh->mesh_dimension(), X_fe_type);
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, d_mesh->mesh_dimension(), THIRD);
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<std::vector<double>>& phi = fe->get_phi();
    const std::vector<double>& JxW = fe->get_JxW();

    // Only changes are needed where the structure lives
    const int level_num = d_fe_data_manager->getFinestPatchLevelNumber();
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_num);

    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        for (const auto& elem_idx_pair : d_elem_idx_map_vec[level_num])
        {
            const Pointer<Patch<NDIM>>& patch = level->getPatch(elem_idx_pair.first.d_patch_num);
            const CellIndex<NDIM>& idx = elem_idx_pair.first.d_idx;
            const std::vector<std::pair<Elem*, std::unique_ptr<Elem>>>& elem_pairs = elem_idx_pair.second;
            const std::vector<std::vector<libMesh::Point>>& intersection_point_vec =
                d_intersection_idx_map_vec[level_num][PatchIndexPair(patch, idx)];

            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(d_area_idx);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
            Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);

            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

            for (unsigned int i = 0; i < elem_pairs.size(); ++i)
            {
                const Elem* const old_elem = elem_pairs[i].first;
                const Elem* const new_elem = elem_pairs[i].second.get();
                const std::vector<libMesh::Point>& intersection_points = intersection_point_vec[i];

                std::vector<dof_id_type> fl_dofs, sf_dofs, Q_dofs;
                boost::multi_array<double, 2> x_node;
                boost::multi_array<double, 1> Q_node;
                std::vector<boost::multi_array<double, 1>> fl_node(d_fl_names.size()), sf_node(d_sf_names.size());
                std::vector<double> sf_vals(d_sf_names.size()), fl_vals(d_fl_names.size());

                const auto& X_dof_indices = X_dof_map_cache.dof_indices(old_elem);
                IBTK::get_values_for_interpolation(x_node, *X_petsc_vec, X_local_soln, X_dof_indices);

                Q_dof_map.dof_indices(old_elem, Q_dofs);
                IBTK::get_values_for_interpolation(Q_node, *Q_vec, Q_dofs);
                for (unsigned int l = 0; l < d_fl_names.size(); ++l)
                {
                    fl_dof_maps[l]->dof_indices(old_elem, fl_dofs);
                    IBTK::get_values_for_interpolation(fl_node[l], *fl_vecs[l], fl_dofs);
                }
                for (unsigned int l = 0; l < d_sf_names.size(); ++l)
                {
                    sf_dof_maps[l]->dof_indices(old_elem, sf_dofs);
                    IBTK::get_values_for_interpolation(sf_node[l], *sf_vecs[l], sf_dofs);
                }
                // We need to interpolate our solution to the new element's nodes
                std::array<double, 2> Q_soln_on_new_elem;
                std::vector<std::array<double, 2>> sf_soln_on_new_elem(d_sf_names.size()),
                    fl_soln_on_new_elem(d_fl_names.size());
                fe->reinit(old_elem, &intersection_points);
                for (unsigned int l = 0; l < 2; ++l)
                {
                    Q_soln_on_new_elem[l] = IBTK::interpolate(l, Q_node, phi);
                    for (unsigned int k = 0; k < d_sf_names.size(); ++k)
                        sf_soln_on_new_elem[k][l] = IBTK::interpolate(l, sf_node[k], phi);
                    for (unsigned int k = 0; k < d_fl_names.size(); ++k)
                        fl_soln_on_new_elem[k][l] = IBTK::interpolate(l, fl_node[k], phi);
                }
                // Then we need to integrate
                fe->reinit(new_elem);
                double a = 0.0, g = 0.0;
                double area = 0.0;
                for (unsigned int qp = 0; qp < JxW.size(); ++qp)
                {
                    double Q_val = 0.0;
                    std::fill(sf_vals.begin(), sf_vals.end(), 0.0);
                    std::fill(fl_vals.begin(), fl_vals.end(), 0.0);
                    for (int n = 0; n < 2; ++n)
                    {
                        Q_val += Q_soln_on_new_elem[n] * phi[n][qp];
                        for (unsigned int l = 0; l < d_fl_names.size(); ++l)
                            fl_vals[l] += fl_soln_on_new_elem[l][n] * phi[n][qp];
                        for (unsigned int l = 0; l < d_sf_names.size(); ++l)
                            sf_vals[l] += sf_soln_on_new_elem[l][n] * phi[n][qp];
                    }
                    a += d_a_fcn(Q_val, fl_vals, sf_vals, time, d_fcn_ctx) * JxW[qp];
                    area += JxW[qp];
                    if (!d_homogeneous_bdry) g += d_g_fcn(Q_val, fl_vals, sf_vals, time, d_fcn_ctx) * JxW[qp];
                }

                double cell_volume = dx[0] * dx[1] * (*vol_data)(idx);
                if (cell_volume <= 0.0)
                {
                    plog << "Found intersection with zero cell volume.\n";
                    plog << "On index: " << idx << "\n";
                    plog << "Ignoring contribution.\n";
                    continue;
                }
                if (!d_homogeneous_bdry) (*R_data)(idx) += pre_fac * g / cell_volume;
                (*R_data)(idx) -= pre_fac * a / cell_volume;
            }
        }
    }
    X_petsc_vec->restore_array();
    LS_TIMER_STOP(t_applyBoundaryCondition);
}

void
SBBoundaryConditions::interpolateToBoundary(const int Q_idx,
                                            const std::string& Q_sys_name,
                                            Pointer<PatchHierarchy<NDIM>> hierarchy,
                                            const double current_time)
{
    LS_TIMER_START(t_interpolateToBoundary);
    // First ensure we've filled ghost cells
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp(1);
    ghost_cell_comp[0] =
        ITC(d_scr_idx, Q_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr);
    HierarchyGhostCellInterpolation hier_ghost_cell;
    hier_ghost_cell.initializeOperatorState(ghost_cell_comp, hierarchy);
    hier_ghost_cell.fillData(current_time);

    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();
    System& Q_system = eq_sys->get_system(Q_sys_name);
    const unsigned int n_vars = Q_system.n_vars();
    const DofMap& Q_dof_map = Q_system.get_dof_map();
    System& X_system = eq_sys->get_system(d_fe_data_manager->COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();

    NumericVector<double>* X_vec = X_system.current_local_solution.get();
    NumericVector<double>* Q_vec = Q_system.solution.get();

    auto X_petsc_vec = dynamic_cast<PetscVector<double>*>(X_vec);
    TBOX_ASSERT(X_petsc_vec != nullptr);
    const double* const X_local_soln = X_petsc_vec->get_array_read();

    Q_vec->zero();

    // Loop over patches and interpolate solution to the boundary
    // Assume we are only doing this on the finest level
    // TODO: This stuff only changes in between timesteps. We should cache this information.
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(hierarchy->getFinestLevelNumber());
    const std::vector<std::vector<Node*>>& active_patch_node_map = d_fe_data_manager->getActivePatchNodeMap();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const int local_patch_num = patch->getPatchNumber();

        const std::vector<Node*>& patch_nodes = active_patch_node_map[local_patch_num];
        const size_t num_active_patch_nodes = patch_nodes.size();
        if (!num_active_patch_nodes) continue;

        const Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const patch_x_low = pgeom->getXLower();
        const double* const patch_x_up = pgeom->getXUpper();
        std::array<bool, NDIM> touches_upper_regular_bdry;
        for (int d = 0; d < NDIM; ++d) touches_upper_regular_bdry[d] = pgeom->getTouchesRegularBoundary(d, 1);

        // Store the value of X at the nodes that are inside the current patch
        std::vector<dof_id_type> Q_node_idxs;
        std::vector<double> Q_node, X_node;
        Q_node_idxs.reserve(n_vars * num_active_patch_nodes);
        Q_node.reserve(n_vars * num_active_patch_nodes);
        X_node.reserve(NDIM * num_active_patch_nodes);
        std::vector<dof_id_type> Q_idxs, X_idxs;
        IBTK::Point X;
        for (unsigned int k = 0; k < num_active_patch_nodes; ++k)
        {
            const Node* const n = patch_nodes[k];
            bool inside_patch = true;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                IBTK::get_nodal_dof_indices(X_dof_map, n, d, X_idxs);
                X[d] = X_local_soln[X_petsc_vec->map_global_to_local_index(X_idxs[0])];
                inside_patch = inside_patch && (X[d] >= patch_x_low[d]) &&
                               ((X[d] < patch_x_up[d]) || (touches_upper_regular_bdry[d] && X[d] <= patch_x_up[d]));
            }
            if (inside_patch)
            {
                Q_node.resize(Q_node.size() + n_vars, 0.0);
                X_node.insert(X_node.end(), &X[0], &X[0] + NDIM);
                for (unsigned int i = 0; i < n_vars; ++i)
                {
                    IBTK::get_nodal_dof_indices(Q_dof_map, n, i, Q_idxs);
                    Q_node_idxs.insert(Q_node_idxs.end(), Q_idxs.begin(), Q_idxs.end());
                }
            }
        }

        TBOX_ASSERT(Q_node.size() <= n_vars * num_active_patch_nodes);
        TBOX_ASSERT(X_node.size() <= NDIM * num_active_patch_nodes);
        TBOX_ASSERT(Q_node_idxs.size() <= n_vars * num_active_patch_nodes);

        if (Q_node.empty()) continue;

        // Now we can interpolate from the fluid to the structure.
        Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(d_scr_idx);
        for (size_t i = 0; i < Q_node.size(); ++i)
        {
            // Use a MLS linear approximation to evaluate data on structure
            const CellIndex<NDIM> cell_idx = IndexUtilities::getCellIndex(&X_node[NDIM * i], pgeom, patch->getBox());
            const CellIndex<NDIM>& idx_low = patch->getBox().lower();
            VectorNd x_loc = {
                X_node[NDIM * i],
                X_node[NDIM * i + 1]
#if (NDIM == 3)
                ,
                X_node[NDIM * i + 2]
#endif
            };
            Q_node[i] = d_rbf_reconstruct.reconstructOnIndex(x_loc, cell_idx, *Q_data, patch);
        }
        Q_vec->add_vector(Q_node, Q_node_idxs);
    }
    X_petsc_vec->restore_array();
    Q_vec->close();
    LS_TIMER_STOP(t_interpolateToBoundary);
}

bool
SBBoundaryConditions::findIntersection(libMesh::Point& p, Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q)
{
    bool found_intersection = false;
    switch (elem->type())
    {
    case libMesh::EDGE2:
    {
        // Use linear interpolation
        // Plane through r in q direction:
        // p = r + t * q
        // Plane through two element points p0, p1
        // p = 0.5*(1+u)*p0 + 0.5*(1-u)*p1
        // Set equal and solve for u and t.
        // Note that since q is aligned with a grid axis, we can solve for u first, then find t later
        // Solve for u via a * u + b = 0
        // with a = 0.5 * (p0 - p1)
        //      b = 0.5 * (p0 + p1) - r
        const libMesh::Point& p0 = elem->point(0);
        const libMesh::Point& p1 = elem->point(1);
        const int search_dir = q(0) == 0.0 ? 1 : 0;
        const int trans_dir = (search_dir + 1) % NDIM;
        double a = 0.5 * (p0(trans_dir) - p1(trans_dir));
        double b = 0.5 * (p0(trans_dir) + p1(trans_dir)) - r(trans_dir);
        const double u = -b / a;
        // Determine if this intersection is on the interior of the element
        // This means that u is between -1 and 1
        if (u >= -1.0 && u <= 1.0)
        {
            // Now determine if intersection occurs on axis
            // This amounts to t being between -0.5 and 0.5
            double p_search = 0.5 * p0(search_dir) * (1.0 + u) + 0.5 * (1.0 - u) * p1(search_dir);
            double t = (p_search - r(search_dir)) / q(search_dir);
            if (t >= -0.5 && t <= 0.5)
            {
                // We've found an intersection on this axis
                p = 0.5 * (1.0 + u) * p0 + 0.5 * (1.0 - u) * p1;
                found_intersection = true;
            }
        }
        break;
    }
    default:
        TBOX_ERROR("Unknown element.\n");
    }
    return found_intersection;
}

void
SBBoundaryConditions::clearCache()
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        d_node_idx_map_vec[ln].clear();
        d_elem_idx_map_vec[ln].clear();
        d_intersection_idx_map_vec[ln].clear();
    }
}

void
SBBoundaryConditions::cacheBoundaryData()
{
    LS_TIMER_START(t_cache_boundary_data);
    d_node_idx_map_vec.resize(d_hierarchy->getFinestLevelNumber() + 1);
    d_elem_idx_map_vec.resize(d_hierarchy->getFinestLevelNumber() + 1);
    d_intersection_idx_map_vec.resize(d_hierarchy->getFinestLevelNumber() + 1);
    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();

    System& X_system = eq_sys->get_system(d_fe_data_manager->COORDINATES_SYSTEM_NAME);
    DofMap& X_dof_map = X_system.get_dof_map();
    FEType X_fe_type = X_dof_map.variable_type(0);
    NumericVector<double>* X_vec = X_system.solution.get();
    auto X_petsc_vec = dynamic_cast<PetscVector<double>*>(X_vec);
    TBOX_ASSERT(X_petsc_vec != nullptr);
    const double* const X_local_soln = X_petsc_vec->get_array_read();
    FEDataManager::SystemDofMapCache& X_dof_map_cache =
        *d_fe_data_manager->getDofMapCache(d_fe_data_manager->COORDINATES_SYSTEM_NAME);

    // Only changes are needed where the structure lives
    const int level_num = d_fe_data_manager->getFinestPatchLevelNumber();
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);
    const Pointer<CartesianGridGeometry<NDIM>> grid_geom = level->getGridGeometry();
    VectorValue<double> n;
    IBTK::Point x_min, x_max;
    const std::vector<std::vector<Elem*>>& active_patch_element_map = d_fe_data_manager->getActivePatchElementMap();

    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const std::vector<Elem*>& patch_elems = active_patch_element_map[patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(d_area_idx);
        Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
        Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);

        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

        std::vector<dof_id_type> fl_dofs, sf_dofs, Q_dofs;
        boost::multi_array<double, 2> x_node;
        boost::multi_array<double, 1> Q_node;
        std::vector<boost::multi_array<double, 1>> fl_node(d_fl_names.size()), sf_node(d_sf_names.size());
        std::vector<double> sf_vals(d_sf_names.size()), fl_vals(d_fl_names.size());
        for (const auto& elem : patch_elems)
        {
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            IBTK::get_values_for_interpolation(x_node, *X_petsc_vec, X_local_soln, X_dof_indices);

            // Simple check if element is completely within grid cell.
            const unsigned int n_node = elem->n_nodes();
            std::vector<hier::Index<NDIM>> elem_idx_nodes(n_node);
            for (unsigned int k = 0; k < n_node; ++k)
            {
                const Node& node = elem->node_ref(k);
                elem_idx_nodes[k] = IndexUtilities::getCellIndex(&node(0), grid_geom, level->getRatio());
            }
            // Check if all indices are the same
            if (std::adjacent_find(elem_idx_nodes.begin(),
                                   elem_idx_nodes.end(),
                                   std::not_equal_to<hier::Index<NDIM>>()) == elem_idx_nodes.end())
            {
                // Element is entirely contained in cell.
                // Store element and continue to next element
                PatchIndexPair p_idx(patch, CellIndex<NDIM>(elem_idx_nodes[0]));
                // Create copy of element
                std::unique_ptr<Elem> new_elem = Elem::build(elem->type());
                std::array<std::unique_ptr<Node>, 2> new_elem_nodes;
                for (int n = 0; n < elem->n_nodes(); ++n)
                {
                    new_elem_nodes[n] = libmesh_make_unique<Node>(*elem->get_nodes()[n], n);
                    new_elem->set_id(0);
                    new_elem->set_node(n) = new_elem_nodes[n].get();
                    d_node_idx_map_vec[level_num][p_idx].push_back(std::move(new_elem_nodes[n]));
                }
                d_elem_idx_map_vec[level_num][p_idx].push_back(std::make_pair(elem, std::move(new_elem)));
                continue;
            }

            std::vector<libMesh::Point> X_node_cache(n_node), x_node_cache(n_node);
            x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
            x_max = IBTK::Point::Constant(-std::numeric_limits<double>::max());
            for (unsigned int k = 0; k < n_node; ++k)
            {
                X_node_cache[k] = elem->point(k);
                libMesh::Point& x = x_node_cache[k];
                for (unsigned int d = 0; d < NDIM; ++d) x(d) = x_node[k][d];
                // Perturb the mesh so that we keep FE nodes away from cell edges / nodes
                // Therefore we don't have to worry about nodes being on a cell edge
                if (d_perturb_nodes)
                {
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int i_s = std::floor((x(d) - x_lower[d]) / dx[d]) + patch_lower[d];
                        for (int shift = 0; shift <= 2; ++shift)
                        {
                            const double x_s =
                                x_lower[d] + dx[d] * (static_cast<double>(i_s - patch_lower[d]) + 0.5 * shift);
                            const double tol = 1.0e-8 * dx[d];
                            if (x(d) <= x_s) x(d) = std::min(x_s - tol, x(d));
                            if (x(d) >= x_s) x(d) = std::max(x_s + tol, x(d));
                        }
                    }
                }

                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x_min[d] = std::min(x_min[d], x(d));
                    x_max[d] = std::max(x_max[d], x(d));
                }
                elem->point(k) = x;
            }
            Box<NDIM> box(IndexUtilities::getCellIndex(&x_min[0], grid_geom, level->getRatio()),
                          IndexUtilities::getCellIndex(&x_max[0], grid_geom, level->getRatio()));
            box.grow(1);
            box = box * patch->getBox();

            // We have the bounding box of the element. Now loop over coordinate directions and look for intersections
            // with the background grid.
            for (BoxIterator<NDIM> b(box); b; b++)
            {
                const hier::Index<NDIM>& i_c = b();
                // We have the index of the box. Each box should have zero or two intersections
                std::vector<libMesh::Point> intersection_points(0);
                for (int upper_lower = 0; upper_lower < 2; ++upper_lower)
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        VectorValue<double> q;
#if (NDIM == 2)
                        q((axis + 1) % NDIM) = dx[(axis + 1) % NDIM];
#endif
                        libMesh::Point r;
                        for (int d = 0; d < NDIM; ++d)
                            r(d) = x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) +
                                                         (d == axis ? (upper_lower == 1 ? 1.0 : 0.0) : 0.5));

                        libMesh::Point p;
                        // An element may intersect zero or one times with a cell edge.
                        if (findIntersection(p, elem, r, q)) intersection_points.push_back(p);
                    }
                }

                // An element may have zero, one, or two intersections with a cell.
                // Note we've already accounted for when an element is contained within a cell.
                if (intersection_points.size() == 0) continue;
                TBOX_ASSERT(intersection_points.size() <= 2);

                // Create a new element
                std::unique_ptr<Elem> new_elem = Elem::build(EDGE2);
                if (intersection_points.size() == 1)
                {
                    for (unsigned int k = 0; k < n_node; ++k)
                    {
                        libMesh::Point xn;
                        for (int d = 0; d < NDIM; ++d) xn(d) = elem->point(k)(d);
                        const hier::Index<NDIM>& n_idx =
                            IndexUtilities::getCellIndex(&xn(0), grid_geom, level->getRatio());
                        if (n_idx == i_c)
                        {
                            // Check if we already have this point accounted for. Note this can happen when a node is
                            // EXACTLY on a cell face or node.
                            if (intersection_points[0] == xn) continue;
                            intersection_points.push_back(xn);
                            break;
                        }
                    }
                }
                // At this point, if we still only have one intersection point, our node is on a face, and we can skip
                // this index.
                if (intersection_points.size() == 1) continue;
                TBOX_ASSERT(intersection_points.size() == 2);
                PatchIndexPair p_idx(patch, CellIndex<NDIM>(i_c));
                int cur_end = d_node_idx_map_vec[level_num][p_idx].size();
                for (int i = 0; i < 2; ++i)
                {
                    d_node_idx_map_vec[level_num][p_idx].push_back(
                        libmesh_make_unique<Node>(intersection_points[i], i));
                    new_elem->set_id(0);
                    new_elem->set_node(i) = d_node_idx_map_vec[level_num][p_idx][i + cur_end].get();
                }
                // We have new element. Add it to the list.
                d_intersection_idx_map_vec[level_num][p_idx].push_back(std::move(intersection_points));
                d_elem_idx_map_vec[level_num][p_idx].push_back(std::make_pair(elem, std::move(new_elem)));
            }
        }
    }
    LS_TIMER_STOP(t_cache_boundary_data);
}
} // namespace LS
