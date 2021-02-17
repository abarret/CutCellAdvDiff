#include "ibtk/IndexUtilities.h"

#include "LS/SBIntegrator.h"
#include "LS/utility_functions.h"

#include "libmesh/elem_cutter.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/transient_system.h"

#include <boost/multi_array.hpp>

namespace
{
static Timer* t_integrateHierarchy = nullptr;
static Timer* t_interpolateToBoundary = nullptr;
} // namespace

namespace LS
{
SBIntegrator::SBIntegrator(std::string object_name,
                           Pointer<Database> input_db,
                           Mesh* mesh,
                           FEDataManager* fe_data_manager)
    : d_object_name(std::move(object_name)),
      d_mesh(mesh),
      d_fe_data_manager(fe_data_manager),
      d_scr_var(new CellVariable<NDIM, double>(d_object_name + "::SCR"))
{
    d_perturb_nodes = input_db->getBool("perturb_nodes");
    d_rbf_reconstruct.setStencilWidth(input_db->getInteger("stencil_width"));

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_scr_idx = var_db->registerVariableAndContext(
        d_scr_var, var_db->getContext(d_object_name + "::SCR"), d_rbf_reconstruct.getStencilWidth());

    IBTK_DO_ONCE(t_integrateHierarchy = TimerManager::getManager()->getTimer("LS::SBIntegrator::integrateHierarchy()");
                 t_interpolateToBoundary =
                     TimerManager::getManager()->getTimer("LS::SBIntegrator::interpolateToBoundary()"););
}

void
SBIntegrator::registerFluidConcentration(Pointer<CellVariable<NDIM, double>> fl_var)
{
    if (d_fe_eqs_initialized)
        TBOX_ERROR(d_object_name + ": can't register a fluid variable after equation systems have been initialized.");
    TBOX_ASSERT(fl_var);
    if (std::find(d_fl_vars.begin(), d_fl_vars.end(), fl_var) == d_fl_vars.end())
    {
        d_fl_vars.push_back(fl_var);
        d_fl_names.push_back(fl_var->getName());
    }
}

void
SBIntegrator::registerFluidConcentration(const std::vector<Pointer<CellVariable<NDIM, double>>>& fl_vars)
{
    for (const auto& fl_var : fl_vars) registerFluidConcentration(fl_var);
}

void
SBIntegrator::registerSurfaceConcentration(std::string surface_name)
{
    if (d_fe_eqs_initialized)
        TBOX_ERROR(d_object_name + ": can't register a surface variable after equation systems have been initialized.");
    if (std::find(d_sf_names.begin(), d_sf_names.end(), surface_name) == d_sf_names.end())
        d_sf_names.push_back(std::move(surface_name));
}

void
SBIntegrator::registerSurfaceConcentration(const std::vector<std::string>& surface_names)
{
    for (const auto& surface_name : surface_names) registerSurfaceConcentration(surface_name);
}

void
SBIntegrator::registerFluidSurfaceDependence(const std::string& surface_name,
                                             Pointer<CellVariable<NDIM, double>> fl_var)
{
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), surface_name) != d_sf_names.end());
    const auto& fl_it = std::find(d_fl_vars.begin(), d_fl_vars.end(), fl_var);
    TBOX_ASSERT(fl_it != d_fl_vars.end());
    const unsigned int l = std::distance(d_fl_vars.begin(), fl_it);
    const std::string& fl_name = d_fl_names[l];
    std::vector<std::string>& fl_vars_vec = d_sf_fl_map[surface_name];
    if (std::find(fl_vars_vec.begin(), fl_vars_vec.end(), fl_name) == fl_vars_vec.end()) fl_vars_vec.push_back(fl_name);
}

void
SBIntegrator::registerSurfaceSurfaceDependence(const std::string& part1_name, const std::string& part2_name)
{
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), part1_name) != d_sf_names.end());
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), part2_name) != d_sf_names.end());
    std::vector<std::string>& sf_names_vec = d_sf_sf_map[part1_name];
    if (std::find(sf_names_vec.begin(), sf_names_vec.end(), part2_name) != sf_names_vec.end())
        sf_names_vec.push_back(part2_name);
}

void
SBIntegrator::registerSurfaceReactionFunction(const std::string& surface_name,
                                              ReactionFcn fcn,
                                              void* ctx /* = nullptr */)
{
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), surface_name) != d_sf_names.end());
    d_sf_reaction_fcn_map[surface_name] = fcn;
    d_sf_ctx_map[surface_name] = ctx;
}

void
SBIntegrator::initializeFEEquationSystems()
{
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();

    if (from_restart)
    {
        TBOX_ERROR("Restart not currently supported!\n\n");
    }
    else
    {
        for (const auto& sf_name : d_sf_names)
        {
            auto& surface_sys = eq_sys->add_system<TransientExplicitSystem>(sf_name);
            surface_sys.add_variable(sf_name, FEType());
            auto& surface_sys_scr = eq_sys->add_system<ExplicitSystem>(sf_name + "_scratch");
            surface_sys_scr.add_variable(sf_name + "_scratch", FEType());
        }

        for (const auto& fl_name : d_fl_names)
        {
            auto& fluid_sys = eq_sys->add_system<ExplicitSystem>(fl_name);
            fluid_sys.add_variable(fl_name, FEType());
        }
    }
}

void
SBIntegrator::setLSData(const int ls_idx, const int vol_idx, Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    d_ls_idx = ls_idx;
    d_vol_idx = vol_idx;
    d_hierarchy = hierarchy;
    d_rbf_reconstruct.setLSData(ls_idx, vol_idx);
    d_rbf_reconstruct.setPatchHierarchy(hierarchy);
}

void
SBIntegrator::integrateHierarchy(Pointer<VariableContext> ctx,
                                 const double current_time,
                                 const double new_time,
                                 const int cycle_num)
{
    plog << "On cycle num: " << cycle_num << "\n";
    // Integrate hierarchy. First loop through, we compute
    // C^tilde = C^n + dt * f(C^n, c_f^n)
    // C^n+1 = C^n + dt * f(C^tilde, c_f^n)
    // Any extra cycles
    // C^tilde = C^n + dt * f(C^n, c_f^n)
    // C^n+1 = C^n + dt * f(C^tilde, c_f^{n+1/2})
    // We assume that the fluid ctx will gives us c_f^{n+1/2} and we will store C^tilde from previous cycles.
    LS_TIMER_START(t_integrateHierarchy);
    // First we need to interpolate all fluid concentrations to the boundary.
    // TODO: We assume that the fluid variable has appropriate data filled with the given context. Should we check this?
    for (unsigned int l = 0; l < d_fl_names.size(); ++l)
    {
        interpolateToBoundary(d_fl_vars[l], ctx, d_hierarchy);
    }
    // Now we actually solve the ODE.
    for (int cn = 0; cn < 2; ++cn)
    {
        // If the global cycle number is greater than one, then we don't need to do first step again
        if (cycle_num > 0) ++cn;
        for (const auto& sf_name : d_sf_names)
        {
            const double dt = new_time - current_time;
            // Solve ODE on surface
            EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();
            auto& sf_base_sys = eq_sys->get_system<TransientExplicitSystem>(sf_name);
            DofMap& sf_base_dof_map = sf_base_sys.get_dof_map();
            NumericVector<double>* sf_base_cur_vec = sf_base_sys.solution.get();
            NumericVector<double>* sf_base_old_vec = sf_base_sys.old_local_solution.get();

            auto& sf_scr_sys = eq_sys->get_system<ExplicitSystem>(sf_name + "_scratch");
            DofMap& sf_scr_dof_map = sf_scr_sys.get_dof_map();
            NumericVector<double>* sf_base_scr_vec = sf_scr_sys.solution.get();

            // Get the NumericVector for all associated systems.
            // TODO: We should check that the dof_maps are all the same, otherwise we need to grab them too.
            const std::vector<std::string>& sf_sf_vec = d_sf_sf_map[sf_name];
            const std::vector<std::string>& sf_fl_vec = d_sf_fl_map[sf_name];
            std::vector<NumericVector<double>*> sf_cur_vecs, fl_vecs;
            std::vector<DofMap*> sf_dof_maps, fl_dof_maps;
            for (const auto& fl_name : sf_fl_vec)
            {
                System& fl_sys = eq_sys->get_system(fl_name);
                fl_dof_maps.push_back(&fl_sys.get_dof_map());
                fl_vecs.push_back(fl_sys.solution.get());
            }
            for (const auto& sf_name : sf_sf_vec)
            {
                if (cn == 0)
                {
                    auto& sf_sys = eq_sys->get_system<TransientExplicitSystem>(sf_name);
                    sf_dof_maps.push_back(&sf_sys.get_dof_map());
                    sf_cur_vecs.push_back(sf_sys.old_local_solution.get());
                }
                else
                {
                    auto& sf_sys = eq_sys->get_system<ExplicitSystem>(sf_name + "_scratch");
                    sf_dof_maps.push_back(&sf_sys.get_dof_map());
                    sf_cur_vecs.push_back(sf_sys.solution.get());
                }
            }

            // Assume our structure is given on finest level.
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(d_hierarchy->getFinestLevelNumber());
            const std::vector<std::vector<Node*>>& active_patch_node_map = d_fe_data_manager->getActivePatchNodeMap();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const int patch_num = patch->getPatchNumber();
                const std::vector<Node*>& active_patch_nodes = active_patch_node_map[patch_num];

                for (const auto& node : active_patch_nodes)
                {
                    // Integrate solution to new value. Use midpoint-rule
                    std::vector<dof_id_type> fl_dofs, sf_dofs, sf_base_dofs, sf_scr_dofs;
                    std::vector<double> fl_vals, sf_cur_vals;
                    std::vector<double> sf_base_cur_vals, sf_base_new_vals;
                    for (unsigned int l = 0; l < sf_sf_vec.size(); ++l)
                    {
                        IBTK::get_nodal_dof_indices(*sf_dof_maps[l], node, 0, sf_dofs);
                        sf_cur_vals.push_back((*sf_cur_vecs[l])(sf_dofs[0]));
                    }
                    for (unsigned int l = 0; l < sf_fl_vec.size(); ++l)
                    {
                        IBTK::get_nodal_dof_indices(*fl_dof_maps[l], node, 0, fl_dofs);
                        fl_vals.push_back((*fl_vecs[l])(fl_dofs[0]));
                    }
                    IBTK::get_nodal_dof_indices(sf_base_dof_map, node, 0, sf_base_dofs);
                    IBTK::get_nodal_dof_indices(sf_scr_dof_map, node, 0, sf_scr_dofs);
                    const double sf_old_val = (*sf_base_old_vec)(sf_base_dofs[0]);
                    double sf_cur_val = (*sf_base_cur_vec)(sf_base_dofs[0]);
                    double sf_scr_val = (*sf_base_scr_vec)(sf_scr_dofs[0]);
                    if (cn == 0)
                    {
                        // We need to fill scratch values
                        sf_scr_val =
                            sf_old_val + 0.5 * dt *
                                             d_sf_reaction_fcn_map[sf_name](
                                                 sf_old_val, fl_vals, sf_cur_vals, current_time, d_sf_ctx_map[sf_name]);
                        sf_base_scr_vec->set(sf_scr_dofs[0], sf_scr_val);
                    }
                    else
                    {
                        // We need to fill in new values
                        sf_cur_val =
                            sf_old_val +
                            dt * d_sf_reaction_fcn_map[sf_name](
                                     sf_scr_val, fl_vals, sf_cur_vals, current_time + 0.5 * dt, d_sf_ctx_map[sf_name]);
                        sf_base_cur_vec->set(sf_base_dofs[0], sf_cur_val);
                    }
                }
            }
            sf_base_cur_vec->close();
            sf_base_scr_vec->close();
        }
    }
    LS_TIMER_STOP(t_integrateHierarchy);
}

void
SBIntegrator::beginTimestepping(const double /*current_time*/, const double /*new_time*/, bool update_time_vals)
{
    // Overwrite old surface concentration
    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();
    if (update_time_vals)
    {
        for (const auto& sf_name : d_sf_names)
        {
            TransientExplicitSystem& sf_sys = eq_sys->get_system<TransientExplicitSystem>(sf_name);
            *sf_sys.older_local_solution = *sf_sys.old_local_solution;
            *sf_sys.old_local_solution = *sf_sys.current_local_solution;
            sf_sys.update();
        }
    }

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_scr_idx)) level->allocatePatchData(d_scr_idx);
    }
}

void
SBIntegrator::endTimestepping(const double /*current_time*/, const double /*new_time*/, bool /*update_time_vals*/)
{
    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();
    for (const auto& sf_name : d_sf_names)
    {
        TransientExplicitSystem& sf_sys = eq_sys->get_system<TransientExplicitSystem>(sf_name);
        sf_sys.update();
    }

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_scr_idx)) level->deallocatePatchData(d_scr_idx);
    }

    d_rbf_reconstruct.clearCache();
}

void
SBIntegrator::interpolateToBoundary(Pointer<CellVariable<NDIM, double>> fl_var,
                                    Pointer<VariableContext> ctx,
                                    const double time)
{
    LS_TIMER_START(t_interpolateToBoundary);
    const auto& fl_it = std::find(d_fl_vars.begin(), d_fl_vars.end(), fl_var);
    TBOX_ASSERT(fl_it != d_fl_vars.end());
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int fl_idx = var_db->mapVariableAndContextToIndex(fl_var, ctx);
    TBOX_ASSERT(fl_idx != IBTK::invalid_index);
    const int l = std::distance(d_fl_vars.begin(), fl_it);
    const std::string& fl_name = d_fl_names[l];
    // First ensure we've filled ghost cells
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp(1);
    ghost_cell_comp[0] =
        ITC(d_scr_idx, fl_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr);
    HierarchyGhostCellInterpolation hier_ghost_cell;
    hier_ghost_cell.initializeOperatorState(ghost_cell_comp, d_hierarchy);
    hier_ghost_cell.fillData(time);

    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();
    System& fl_system = eq_sys->get_system(fl_name);
    const unsigned int n_vars = fl_system.n_vars();
    const DofMap& fl_dof_map = fl_system.get_dof_map();
    System& X_system = eq_sys->get_system(d_fe_data_manager->COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();

    NumericVector<double>* X_vec = X_system.current_local_solution.get();
    NumericVector<double>* fl_vec = fl_system.solution.get();

    auto X_petsc_vec = dynamic_cast<PetscVector<double>*>(X_vec);
    TBOX_ASSERT(X_petsc_vec != nullptr);
    const double* const X_local_soln = X_petsc_vec->get_array_read();

    fl_vec->zero();

    // Loop over patches and interpolate solution to the boundary
    // Assume we are only doing this on the finest level
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(d_hierarchy->getFinestLevelNumber());
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
        std::vector<dof_id_type> fl_node_idxs;
        std::vector<double> fl_node, X_node;
        fl_node_idxs.reserve(n_vars * num_active_patch_nodes);
        fl_node.reserve(n_vars * num_active_patch_nodes);
        X_node.reserve(NDIM * num_active_patch_nodes);
        std::vector<dof_id_type> fl_idxs, X_idxs;
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
                fl_node.resize(fl_node.size() + n_vars, 0.0);
                X_node.insert(X_node.end(), &X[0], &X[0] + NDIM);
                for (unsigned int i = 0; i < n_vars; ++i)
                {
                    IBTK::get_nodal_dof_indices(fl_dof_map, n, i, fl_idxs);
                    fl_node_idxs.insert(fl_node_idxs.end(), fl_idxs.begin(), fl_idxs.end());
                }
            }
        }

        TBOX_ASSERT(fl_node.size() <= n_vars * num_active_patch_nodes);
        TBOX_ASSERT(X_node.size() <= NDIM * num_active_patch_nodes);
        TBOX_ASSERT(fl_node_idxs.size() <= n_vars * num_active_patch_nodes);

        if (fl_node.empty()) continue;

        // Now we can interpolate from the fluid to the structure.
        Pointer<CellData<NDIM, double>> fl_data = patch->getPatchData(d_scr_idx);
        Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);
        Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
        for (size_t i = 0; i < fl_node.size(); ++i)
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
            fl_node[i] = d_rbf_reconstruct.reconstructOnIndex(x_loc, cell_idx, *fl_data, patch);
        }
        fl_vec->add_vector(fl_node, fl_node_idxs);
    }
    X_petsc_vec->restore_array();
    fl_vec->close();
    LS_TIMER_STOP(t_interpolateToBoundary);
}
} // namespace LS
