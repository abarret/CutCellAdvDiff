#include "LS/utility_functions.h"

#include "SurfaceBoundaryReactions.h"

#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"

#include <boost/multi_array.hpp>

#include "ibamr/app_namespaces.h"

namespace
{
static Timer* s_apply_timer = nullptr;
}

std::string SurfaceBoundaryReactions::s_surface_sys_name = "SURFACE_CONCENTRATION";
std::string SurfaceBoundaryReactions::s_fluid_sys_name = "FLUID_CONCENTRATION";

SurfaceBoundaryReactions::SurfaceBoundaryReactions(const std::string& object_name,
                                                   Pointer<Database> input_db,
                                                   Mesh* mesh,
                                                   FEDataManager* fe_data_manager)
    : LSCutCellBoundaryConditions(object_name), d_mesh(mesh), d_fe_data_manager(fe_data_manager)
{
    d_k_on = input_db->getDouble("k_on");
    d_k_off = input_db->getDouble("k_off");
    d_D_coef = input_db->getDouble("D");
    d_cb_max = input_db->getDouble("cb_max");
    d_perturb_nodes = input_db->getBool("perturb_nodes");
    d_use_rbfs = input_db->getBool("use_rbfs");

    // Set up equation systems for reactions
    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();
    ExplicitSystem& surface_sys = eq_sys->add_system<ExplicitSystem>(s_surface_sys_name);
    ExplicitSystem& fluid_sys = eq_sys->add_system<ExplicitSystem>(s_fluid_sys_name);
    surface_sys.assemble_before_solve = false;
    surface_sys.assemble();
    fluid_sys.assemble_before_solve = false;
    fluid_sys.assemble();

    surface_sys.add_variable("SurfaceConcentration", FEType());
    fluid_sys.add_variable("FluidConcentration", FEType());

    IBAMR_DO_ONCE(s_apply_timer =
                      TimerManager::getManager()->getTimer("LS::SurfaceBoundaryReactions::applyBoundaryCondition"));
}

void
SurfaceBoundaryReactions::applyBoundaryCondition(Pointer<CellVariable<NDIM, double>> Q_var,
                                                 const int Q_idx,
                                                 Pointer<CellVariable<NDIM, double>> R_var,
                                                 const int R_idx,
                                                 Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                 const double time)
{
    LS_TIMER_START(s_apply_timer);
    TBOX_ASSERT(d_ls_var && d_vol_var && d_area_var);
    TBOX_ASSERT(d_ls_idx > 0 && d_vol_idx > 0 && d_area_idx > 0);

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

    System& Q_fl_sys = eq_sys->get_system(s_fluid_sys_name);
    DofMap& Q_fl_dof_map = Q_fl_sys.get_dof_map();
    FEType Q_fl_fe_type = Q_fl_dof_map.variable_type(0);
    NumericVector<double>* Q_fl_vec = Q_fl_sys.solution.get();

    System& Q_st_sys = eq_sys->get_system(s_surface_sys_name);
    DofMap& Q_st_dof_map = Q_st_sys.get_dof_map();
    FEType Q_st_fe_type = Q_st_dof_map.variable_type(0);
    NumericVector<double>* Q_st_vec = Q_st_sys.solution.get();

    TBOX_ASSERT(Q_st_fe_type == X_fe_type && Q_fl_fe_type == X_fe_type);

    std::unique_ptr<FEBase> fe = FEBase::build(d_mesh->mesh_dimension(), X_fe_type);
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, d_mesh->mesh_dimension(), THIRD);
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<std::vector<double>>& phi = fe->get_phi();
    std::array<const std::vector<std::vector<double>>*, NDIM - 1> dphi_dxi;
    dphi_dxi[0] = &fe->get_dphidxi();
    if (NDIM > 2) dphi_dxi[1] = &fe->get_dphideta();

    // Only changes are needed where the structure lives
    const int level_num = d_fe_data_manager->getFinestPatchLevelNumber();
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_num);
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

        Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
        Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
        Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(d_area_idx);
        Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);

        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

        std::vector<dof_id_type> Q_st_dof_indices;
        boost::multi_array<double, 2> x_node;
        boost::multi_array<double, 1> q_st_node;
        for (const auto& elem : patch_elems)
        {
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            Q_st_dof_map.dof_indices(elem, Q_st_dof_indices);
            IBTK::get_values_for_interpolation(x_node, *X_petsc_vec, X_local_soln, X_dof_indices);
            IBTK::get_values_for_interpolation(q_st_node, *Q_st_vec, Q_st_dof_indices);

            const unsigned int n_node = elem->n_nodes();
            std::vector<libMesh::Point> X_node_cache(n_node), x_node_cache(n_node);
            x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
            x_max = IBTK::Point::Constant(std::numeric_limits<double>::min());
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
                            const double tol = 1.0e-4 * dx[d];
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
                        for (unsigned int d = 0; d < NDIM; ++d)
                            r(d) = x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) +
                                                         (d == axis ? (upper_lower == 1 ? 1.0 : 0.0) : 0.5));

                        libMesh::Point p;
                        // An element may intersect zero or one times with a cell edge.
                        if (findIntersection(p, elem, r, q))
                        {
                            intersection_points.push_back(p);
                        }
                    }
                }

                // An element may have zero, one, or two intersections with a cell.
                TBOX_ASSERT(intersection_points.size() <= 2);
                if (intersection_points.size() == 2)
                {
                    // An element completely pierces this grid cell.
                    // We need to integrate from intersection to intersection.
                    fe->reinit(elem, &intersection_points);
                    double a = 0;
                    double g = 0;
                    for (unsigned int l = 0; l < 2; ++l)
                    {
                        double q_sf_val = IBTK::interpolate(l, q_st_node, phi);
                        a += d_k_on * (d_cb_max - q_sf_val);
                        g -= d_k_off * q_sf_val;
                    }
                    double area = (*area_data)(i_c);
                    double cell_volume = dx[0] * dx[1] * (*vol_data)(i_c);
                    if (cell_volume <= 0.0)
                    {
                        plog << "Found intersection with zero cell volume.\n";
                        plog << "On index: " << i_c << " with intersection points:\n";
                        plog << "P0: " << intersection_points[0] << " and P1: " << intersection_points[1] << "\n";
                        continue;
                    }
                    (*R_data)(i_c) += pre_fac * g * area / cell_volume;
                    (*R_data)(i_c) -= pre_fac * a * (*Q_data)(i_c)*area / cell_volume;
                }
                else if (intersection_points.size() == 1)
                {
                    // An element has one node interior to the cell.
                    // We need to integrate from intersection to node.
                    // First determine the node that is interior to the cell.
                    libMesh::Point x;
                    for (unsigned int k = 0; k < n_node; ++k)
                    {
                        libMesh::Point xn;
                        for (int d = 0; d < NDIM; ++d) xn(d) = x_node[k][d];
                        const Index<NDIM>& n_idx = IndexUtilities::getCellIndex(&xn(0), grid_geom, level->getRatio());
                        if (n_idx == i_c)
                        {
                            x = xn;
                            break;
                        }
                    }
                    intersection_points.push_back(x);
                    fe->reinit(elem, &intersection_points);
                    double a = 0;
                    double g = 0;
                    for (unsigned int l = 0; l < 2; ++l)
                    {
                        double q_sf_val = IBTK::interpolate(l, q_st_node, phi);
                        a += d_k_on * (d_cb_max - q_sf_val);
                        g -= d_k_off * q_sf_val;
                    }
                    double area = (*area_data)(i_c);
                    double cell_volume = dx[0] * dx[1] * (*vol_data)(i_c);
                    if (cell_volume <= 0.0)
                    {
                        plog << "Found intersection with zero cell volume.\n";
                        plog << "On index: " << i_c << " with intersection points:\n";
                        plog << "P0: " << intersection_points[0] << " and P1: " << intersection_points[1] << "\n";
                        continue;
                    }
                    (*R_data)(i_c) += pre_fac * g * area / cell_volume;
                    (*R_data)(i_c) -= pre_fac * a * (*Q_data)(i_c)*area / cell_volume;
                }
            }

            // Restore the element coordinates.
            for (unsigned int k = 0; k < n_node; ++k)
            {
                elem->point(k) = X_node_cache[k];
            }
        }
    }
    Q_st_vec->close();
    X_petsc_vec->restore_array();
    LS_TIMER_STOP(s_apply_timer);
}

void
SurfaceBoundaryReactions::updateSurfaceConcentration(const int fl_idx,
                                                     const double current_time,
                                                     const double new_time,
                                                     Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    interpolateToBoundary(fl_idx, s_fluid_sys_name, hierarchy, current_time);
    const double dt = new_time - current_time;
    auto f = [](double c_fl, double c_sf, double c_sf_max, double k_on, double k_off) -> double {
        return k_on * (c_sf_max - c_sf) * c_fl - k_off * c_sf;
    };
    // Solve ODE on surface
    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();
    System& q_fl_system = eq_sys->get_system(s_fluid_sys_name);
    DofMap& q_fl_dof_map = q_fl_system.get_dof_map();
    System& q_sf_system = eq_sys->get_system(s_surface_sys_name);
    DofMap& q_sf_dof_map = q_sf_system.get_dof_map();

    // Get underlying solution data
    NumericVector<double>* q_fl_vec = q_fl_system.solution.get();
    auto q_fl_petsc_vec = static_cast<PetscVector<double>*>(q_fl_vec);
    const double* const q_fl_local_soln = q_fl_petsc_vec->get_array_read();

    NumericVector<double>* q_sf_vec = q_sf_system.solution.get();

    // Assume we are on the finest level.
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(hierarchy->getFinestLevelNumber());
    const std::vector<std::vector<Node*>>& active_patch_node_map = d_fe_data_manager->getActivePatchNodeMap();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const std::vector<Node*>& active_patch_nodes = active_patch_node_map[patch_num];
        std::vector<dof_id_type> q_fl_dofs, q_sf_dofs;
        std::vector<double> q_sf_vals;

        for (const auto& node : active_patch_nodes)
        {
            // Integrate solution to new value. Use RK2
            IBTK::get_nodal_dof_indices(q_fl_dof_map, node, 0, q_fl_dofs);
            IBTK::get_nodal_dof_indices(q_sf_dof_map, node, 0, q_sf_dofs);
            dof_id_type q_fl_local_dof = q_fl_petsc_vec->map_global_to_local_index(q_fl_dofs[0]);
            q_sf_vec->get(q_sf_dofs, q_sf_vals);
            double q_sf_star =
                q_sf_vals[0] + 0.5 * dt * f(q_fl_local_soln[q_fl_local_dof], q_sf_vals[0], d_cb_max, d_k_on, d_k_off);
            q_sf_vals[0] += dt * f(q_fl_local_soln[q_fl_local_dof], q_sf_star, d_cb_max, d_k_on, d_k_off);
            q_sf_vec->set(q_sf_dofs[0], q_sf_vals[0]);
        }
    }
    q_sf_vec->close();
    q_fl_petsc_vec->restore_array();
}

void
SurfaceBoundaryReactions::spreadToFluid(Pointer<PatchHierarchy<NDIM>> hierarchy)
{
}

void
SurfaceBoundaryReactions::interpolateToBoundary(const int Q_idx,
                                                const std::string& Q_sys_name,
                                                Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                const double current_time)
{
    // First ensure we've filled ghost cells
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp(1);
    ghost_cell_comp[0] = ITC(Q_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR");
    HierarchyGhostCellInterpolation hier_ghost_cell;
    hier_ghost_cell.initializeOperatorState(ghost_cell_comp, hierarchy);
    hier_ghost_cell.fillData(current_time);

    EquationSystems* eq_sys = d_fe_data_manager->getEquationSystems();
    System& Q_system = eq_sys->get_system(Q_sys_name);
    const unsigned int n_vars = Q_system.n_vars();
    const DofMap& Q_dof_map = Q_system.get_dof_map();
    System& X_system = eq_sys->get_system(d_fe_data_manager->COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEType Q_fe_type = Q_dof_map.variable_type(0);
    Order Q_order = Q_dof_map.variable_order(0);

    NumericVector<double>* X_vec = X_system.current_local_solution.get();
    NumericVector<double>* Q_vec = Q_system.solution.get();

    auto X_petsc_vec = dynamic_cast<PetscVector<double>*>(X_vec);
    TBOX_ASSERT(X_petsc_vec != nullptr);
    const double* const X_local_soln = X_petsc_vec->get_array_read();

    Q_vec->zero();

    // Loop over patches and interpolate solution to the boundary
    // Assume we are only doing this on the finest level
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
        const double* const dx = pgeom->getDx();
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
        Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
        Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);
        Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
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
            Q_node[i] = d_use_rbfs ? reconstructRBF(x_loc, cell_idx, *ls_data, *vol_data, *Q_data, patch) :
                                     reconstructMLS(x_loc, cell_idx, *ls_data, *vol_data, *Q_data, patch);
        }
        Q_vec->add_vector(Q_node, Q_node_idxs);
    }
    X_petsc_vec->restore_array();
    Q_vec->close();
}

bool
SurfaceBoundaryReactions::findIntersection(libMesh::Point& p,
                                           Elem* elem,
                                           libMesh::Point r,
                                           libMesh::VectorValue<double> q)
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

double
SurfaceBoundaryReactions::reconstructMLS(const VectorNd& x,
                                         const CellIndex<NDIM>& idx_interp,
                                         const NodeData<NDIM, double>& ls_data,
                                         const CellData<NDIM, double>& vol_data,
                                         const CellData<NDIM, double>& Q_data,
                                         Pointer<Patch<NDIM>> patch)
{
    Box<NDIM> interp_box(idx_interp, idx_interp);
    int box_size = 2;
    int interp_size = NDIM + 1;
    interp_box.grow(box_size);
    TBOX_ASSERT(Q_data.getGhostBox().contains(interp_box));
    TBOX_ASSERT(ls_data.getGhostBox().contains(interp_box));
    TBOX_ASSERT(vol_data.getGhostBox().contains(interp_box));
    std::vector<double> Q_vals;
    std::vector<VectorNd> X_vals;

    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const hier::Index<NDIM>& idx_low = patch->getBox().lower();
    const double* const xlow = pgeom->getXLower();

    for (CellIterator<NDIM> ci(interp_box); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();
        if (vol_data(idx) > 0.0)
        {
            VectorNd x_cent_c = LS::find_cell_centroid(idx, ls_data);
            for (int d = 0; d < NDIM; ++d)
                x_cent_c[d] = xlow[d] + dx[d] * (x_cent_c[d] - static_cast<double>(idx_low(d)));
            Q_vals.push_back(Q_data(idx));
            X_vals.push_back(x_cent_c);
        }
    }
    const int m = Q_vals.size();
    MatrixXd A(MatrixXd::Zero(m, interp_size));
    VectorXd U(VectorXd::Zero(m));
    auto w = [](double r) -> double { return std::exp(-r * r); };
    for (size_t ii = 0; ii < Q_vals.size(); ++ii)
    {
        const VectorNd X = X_vals[ii] - x;
        double weight = std::sqrt(w(X.norm()));
        U(ii) = weight * Q_vals[ii];
        A(ii, 2) = weight * X[1];
        A(ii, 1) = weight * X[0];
        A(ii, 0) = weight;
    }
    VectorXd soln = A.fullPivHouseholderQr().solve(U);
    return soln(0);
}

double
SurfaceBoundaryReactions::reconstructRBF(const VectorNd& x_loc,
                                         const CellIndex<NDIM>& idx_interp,
                                         const NodeData<NDIM, double>& ls_data,
                                         const CellData<NDIM, double>& vol_data,
                                         const CellData<NDIM, double>& Q_data,
                                         Pointer<Patch<NDIM>> patch)
{
    Box<NDIM> box(idx_interp, idx_interp);
    box.grow(2);
#ifndef NDEBUG
    TBOX_ASSERT(ls_data.getGhostBox().contains(box));
    TBOX_ASSERT(Q_data.getGhostBox().contains(box));
    TBOX_ASSERT(vol_data.getGhostBox().contains(box));
#endif

    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_low = pgeom->getXLower();

    auto rbf = [](double r) -> double { return r * r * r; };

    const CellIndex<NDIM>& idx_low = patch->getBox().lower();

    std::vector<double> Q_vals;
    std::vector<VectorNd> X_vals;

    for (CellIterator<NDIM> ci(box); ci; ci++)
    {
        const CellIndex<NDIM>& idx_c = ci();
        if (vol_data(idx_c) > 0.0)
        {
            // Use this point to calculate least squares reconstruction.
            // Find cell center
            VectorNd x_cent_c = LS::find_cell_centroid(idx_c, ls_data);
            for (int d = 0; d < NDIM; ++d)
                x_cent_c[d] = x_low[d] + dx[d] * (x_cent_c(d) - static_cast<double>(idx_low(d)));
            Q_vals.push_back(Q_data(idx_c));
            X_vals.push_back(x_cent_c);
        }
    }
    const int m = Q_vals.size();
    MatrixXd A(MatrixXd::Zero(m, m));
    MatrixXd B(MatrixXd::Zero(m, NDIM + 1));
    VectorXd U(VectorXd::Zero(m + NDIM + 1));
    for (size_t i = 0; i < Q_vals.size(); ++i)
    {
        for (size_t j = 0; j < Q_vals.size(); ++j)
        {
            const VectorNd X = X_vals[i] - X_vals[j];
            A(i, j) = rbf(X.norm());
        }
        B(i, 0) = 1.0;
        for (int d = 0; d < NDIM; ++d) B(i, d + 1) = X_vals[i](d);
        U(i) = Q_vals[i];
    }

    MatrixXd final_mat(MatrixXd::Zero(m + NDIM + 1, m + NDIM + 1));
    final_mat.block(0, 0, m, m) = A;
    final_mat.block(0, m, m, NDIM + 1) = B;
    final_mat.block(m, 0, NDIM + 1, m) = B.transpose();

    VectorXd x = final_mat.fullPivHouseholderQr().solve(U);
    double val = 0.0;
    VectorXd rbf_coefs = x.block(0, 0, m, 1);
    VectorXd poly_coefs = x.block(m, 0, NDIM + 1, 1);
    Vector3d poly_vec = { 1.0, x_loc(0), x_loc(1) };
    for (size_t i = 0; i < X_vals.size(); ++i)
    {
        val += rbf_coefs[i] * rbf((X_vals[i] - x_loc).norm());
    }
    val += poly_coefs.dot(poly_vec);
    return val;
}
