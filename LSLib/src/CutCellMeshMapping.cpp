#include "ibtk/IndexUtilities.h"

#include "LS/CutCellMeshMapping.h"

#include "ibamr/app_namespaces.h"

namespace LS
{
CutCellMeshMapping::CutCellMeshMapping(std::string object_name,
                                       Pointer<Database> input_db,
                                       Mesh* mesh,
                                       FEDataManager* fe_data_manager)
    : d_object_name(std::move(object_name)), d_mesh(mesh), d_fe_data_manager(fe_data_manager)
{
    if (input_db) d_perturb_nodes = input_db->getBool("perturb_nodes");
}

void
CutCellMeshMapping::setLSData(const int ls_idx, const int vol_idx, const int area_idx)
{
    d_ls_idx = ls_idx;
    d_vol_idx = vol_idx;
    d_area_idx = area_idx;
}

void
CutCellMeshMapping::initializeObjectState(Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    if (d_is_initialized) deinitializeObjectState();
    d_hierarchy = hierarchy;

    // Reset mappings
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    d_idx_cut_cell_elems_map_vec.resize(finest_ln + 1);

    d_is_initialized = true;
}

void
CutCellMeshMapping::deinitializeObjectState()
{
    d_idx_cut_cell_elems_map_vec.clear();

    d_is_initialized = false;
}

void
CutCellMeshMapping::generateCutCellMappings()
{
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
                d_idx_cut_cell_elems_map_vec[level_num][p_idx].push_back(
                    CutCellElems(elem, { elem->get_nodes()[0], elem->get_nodes()[1] }));
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
                // Create a new element
                d_idx_cut_cell_elems_map_vec[level_num][p_idx].push_back(CutCellElems(elem, intersection_points));
            }
        }
    }
}

bool
CutCellMeshMapping::findIntersection(libMesh::Point& p, Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q)
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
} // namespace LS
