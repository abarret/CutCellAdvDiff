#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"

#include "LS/SBSurfaceFluidCouplingManager.h"

#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/transient_system.h"

namespace
{
static Timer* t_interpolateToBoundary = nullptr;
}

namespace LS
{
SBSurfaceFluidCouplingManager::SBSurfaceFluidCouplingManager(std::string object_name,
                                                             const Pointer<Database>& input_db,
                                                             FEDataManager* fe_data_manager,
                                                             Mesh* mesh)
    : d_object_name(std::move(object_name)),
      d_fe_data_manager(fe_data_manager),
      d_mesh(mesh),
      d_J_sys_name("Jacobian"),
      d_scr_var(new CellVariable<NDIM, double>(d_object_name + "::SCR"))
{
    d_rbf_reconstruct.setStencilWidth(input_db->getInteger("stencil_width"));

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_scr_idx = var_db->registerVariableAndContext(
        d_scr_var, var_db->getContext(d_object_name + "::SCR"), d_rbf_reconstruct.getStencilWidth());

    IBTK_DO_ONCE(t_interpolateToBoundary =
                     TimerManager::getManager()->getTimer("LS::SBDataManager::interpolateToBoundary()"););
}

SBSurfaceFluidCouplingManager::~SBSurfaceFluidCouplingManager()
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_scr_idx)) level->deallocatePatchData(d_scr_idx);
    }
}

void
SBSurfaceFluidCouplingManager::registerFluidConcentration(Pointer<CellVariable<NDIM, double>> fl_var)
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
SBSurfaceFluidCouplingManager::registerFluidConcentration(
    const std::vector<Pointer<CellVariable<NDIM, double>>>& fl_vars)
{
    for (const auto& fl_var : fl_vars) registerFluidConcentration(fl_var);
}

void
SBSurfaceFluidCouplingManager::registerSurfaceConcentration(std::string surface_name)
{
    if (d_fe_eqs_initialized)
        TBOX_ERROR(d_object_name + ": can't register a surface variable after equation systems have been initialized.");
    if (std::find(d_sf_names.begin(), d_sf_names.end(), surface_name) == d_sf_names.end())
        d_sf_names.push_back(std::move(surface_name));
}

void
SBSurfaceFluidCouplingManager::registerSurfaceConcentration(const std::vector<std::string>& surface_names)
{
    for (const auto& surface_name : surface_names) registerSurfaceConcentration(surface_name);
}

void
SBSurfaceFluidCouplingManager::registerFluidSurfaceDependence(const std::string& sf_name,
                                                              Pointer<CellVariable<NDIM, double>> fl_var)
{
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), sf_name) != d_sf_names.end());
    const auto& fl_it = std::find(d_fl_vars.begin(), d_fl_vars.end(), fl_var);
    TBOX_ASSERT(fl_it != d_fl_vars.end());
    const unsigned int l = std::distance(d_fl_vars.begin(), fl_it);
    const std::string& fl_name = d_fl_names[l];
    std::vector<std::string>& fl_vars_vec = d_sf_fl_map[sf_name];
    if (std::find(fl_vars_vec.begin(), fl_vars_vec.end(), fl_name) == fl_vars_vec.end()) fl_vars_vec.push_back(fl_name);
    std::vector<std::string>& sf_vars_vec = d_fl_sf_map[fl_name];
    if (std::find(sf_vars_vec.begin(), sf_vars_vec.end(), sf_name) == sf_vars_vec.end()) sf_vars_vec.push_back(sf_name);
}

void
SBSurfaceFluidCouplingManager::registerSurfaceSurfaceDependence(const std::string& part1_name,
                                                                const std::string& part2_name)
{
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), part1_name) != d_sf_names.end());
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), part2_name) != d_sf_names.end());
    std::vector<std::string>& sf_names_vec = d_sf_sf_map[part1_name];
    if (std::find(sf_names_vec.begin(), sf_names_vec.end(), part2_name) != sf_names_vec.end())
        sf_names_vec.push_back(part2_name);
}

void
SBSurfaceFluidCouplingManager::registerSurfaceReactionFunction(const std::string& surface_name,
                                                               ReactionFcn fcn,
                                                               void* ctx /* = nullptr */)
{
    TBOX_ASSERT(std::find(d_sf_names.begin(), d_sf_names.end(), surface_name) != d_sf_names.end());
    d_sf_reaction_fcn_ctx_map[surface_name] = std::make_pair(fcn, ctx);
}

void
SBSurfaceFluidCouplingManager::registerFluidBoundaryCondition(const Pointer<CellVariable<NDIM, double>>& fl_var,
                                                              ReactionFcn a_fcn,
                                                              ReactionFcn g_fcn,
                                                              void* ctx)
{
    TBOX_ASSERT(std::find(d_fl_vars.begin(), d_fl_vars.end(), fl_var) != d_fl_vars.end());
    const std::string& fl_name = fl_var->getName();
    TBOX_ASSERT(std::find(d_fl_names.begin(), d_fl_names.end(), fl_name) != d_fl_names.end());
    d_fl_a_g_fcn_map[fl_name] = std::make_tuple(a_fcn, g_fcn, ctx);
}

void
SBSurfaceFluidCouplingManager::initializeFEEquationSystems()
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
        }

        for (const auto& fl_name : d_fl_names)
        {
            auto& fluid_sys = eq_sys->add_system<ExplicitSystem>(fl_name);
            fluid_sys.add_variable(fl_name, FEType());
        }
    }
}

void
SBSurfaceFluidCouplingManager::setLSData(const int ls_idx, const int vol_idx, Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    d_ls_idx = ls_idx;
    d_vol_idx = vol_idx;
    d_hierarchy = hierarchy;
    d_rbf_reconstruct.setLSData(ls_idx, vol_idx);
    d_rbf_reconstruct.setPatchHierarchy(hierarchy);
}

const std::string&
SBSurfaceFluidCouplingManager::interpolateToBoundary(Pointer<CellVariable<NDIM, double>> fl_var,
                                                     const int fl_idx,
                                                     const double time)
{
    LS_TIMER_START(t_interpolateToBoundary);
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_scr_idx)) level->allocatePatchData(d_scr_idx);
    }
    const auto& fl_it = std::find(d_fl_vars.begin(), d_fl_vars.end(), fl_var);
    TBOX_ASSERT(fl_it != d_fl_vars.end());
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
            const CellIndex<NDIM> cell_idx =
                IBTK::IndexUtilities::getCellIndex(&X_node[NDIM * i], pgeom, patch->getBox());
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
    return fl_name;
}

const std::string&
SBSurfaceFluidCouplingManager::updateJacobian()
{
    TBOX_ERROR("Not implemented currently.\n");
    return d_J_sys_name;
}
} // namespace  LS
