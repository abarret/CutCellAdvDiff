#include "ibtk/DebuggingUtilities.h"
#include "ibtk/IBTK_MPI.h"

#include "LS/utility_functions.h"

#include "LSPipeFlow.h"

#include "ibamr/app_namespaces.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define SIGN_SWEEP_FC IBAMR_FC_FUNC(signsweep2dn, SIGNSWEEP2D)
#endif

extern "C"
{
    void SIGN_SWEEP_FC(double* U,
                       const int& U_gcw,
                       const int& ilower0,
                       const int& iupper0,
                       const int& ilower1,
                       const int& iupper1,
                       const double& large_dist,
                       int& n_updates);
}

namespace
{
static Timer* t_updateVolumeAreaSideLS = nullptr;
static Timer* t_findIntersection = nullptr;
} // namespace

namespace LS
{
double LSPipeFlow::s_large_val = std::numeric_limits<double>::max();

LSPipeFlow::LSPipeFlow(std::string object_name,
                       Pointer<PatchHierarchy<NDIM>> hierarchy,
                       MeshBase* lower_mesh,
                       MeshBase* upper_mesh,
                       FEDataManager* lower_manager,
                       FEDataManager* upper_manager,
                       Pointer<Database> input_db)
    : LSFindCellVolume(std::move(object_name), hierarchy),
      d_lower_mesh(lower_mesh),
      d_upper_mesh(upper_mesh),
      d_lower_manager(lower_manager),
      d_upper_manager(upper_manager),
      d_sgn_var(new CellVariable<NDIM, double>(d_object_name + "SGN"))
{
    d_y_low = input_db->getDouble("y_low");
    d_y_up = input_db->getDouble("y_up");
    d_L = input_db->getDouble("l");
    d_theta = input_db->getDouble("theta");

    IBAMR_DO_ONCE(t_updateVolumeAreaSideLS =
                      TimerManager::getManager()->getTimer("LS::LSPipeFlow::updateVolumeAreaSideLS()");
                  t_findIntersection = TimerManager::getManager()->getTimer("LS::LSPipeFlow::findIntersection()"););
    return;
} // Constructor

void
LSPipeFlow::updateVolumeAreaSideLS(int vol_idx,
                                   Pointer<CellVariable<NDIM, double>> /*vol_var*/,
                                   int area_idx,
                                   Pointer<CellVariable<NDIM, double>> /*area_var*/,
                                   int side_idx,
                                   Pointer<SideVariable<NDIM, double>> /*side_var*/,
                                   int phi_idx,
                                   Pointer<NodeVariable<NDIM, double>> phi_var,
                                   double /*data_time*/,
                                   bool extended_box)
{
    LS_TIMER_START(t_updateVolumeAreaSideLS);
    auto dist_up = [this](VectorNd x_pt, double y_p) -> double {
        VectorNd x_int;
        x_int(0) =
            1.0 / std::tan(d_theta) * (x_pt(0) / std::tan(d_theta) + x_pt(1) - y_p) * sin(d_theta) * sin(d_theta);
        x_int(1) = (x_pt(0) / std::tan(d_theta) + x_pt(1) + y_p / (std::tan(d_theta) * std::tan(d_theta))) *
                   sin(d_theta) * sin(d_theta);
        return (x_pt(1) > x_int(1) ? 1.0 : -1.0) * ((x_pt - x_int).norm());
    };
    auto dist_low = [this](VectorNd x_pt, double y_p) -> double {
        VectorNd x_int;
        x_int(0) =
            1.0 / std::tan(d_theta) * (x_pt(0) / std::tan(d_theta) + x_pt(1) - y_p) * sin(d_theta) * sin(d_theta);
        x_int(1) = (x_pt(0) / std::tan(d_theta) + x_pt(1) + y_p / (std::tan(d_theta) * std::tan(d_theta))) *
                   sin(d_theta) * sin(d_theta);
        return (x_pt(1) > x_int(1) ? -1.0 : 1.0) * ((x_pt - x_int).norm());
    };
    TBOX_ASSERT(phi_idx != IBTK::invalid_index);
    TBOX_ASSERT(phi_var);
    HierarchyNodeDataOpsReal<NDIM, double> hier_nc_data_ops(d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    hier_nc_data_ops.setToScalar(phi_idx, s_large_val, false);
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    hier_cc_data_ops.setToScalar(vol_idx, s_large_val, false);
    hier_cc_data_ops.setToScalar(area_idx, 0.0, true);

    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(d_hierarchy->getFinestLevelNumber());

    std::array<FEDataManager*, 2> fe_data_managers = { d_lower_manager, d_upper_manager };
    std::map<hier::Index<NDIM>, std::vector<libMesh::Point>> index_volume_pts_map;
    for (const auto& fe_data_manager : fe_data_managers)
    {
        const std::vector<std::vector<Elem*>>& active_patch_elem_map = fe_data_manager->getActivePatchElementMap();
        const Pointer<CartesianGridGeometry<NDIM>> grid_geom = level->getGridGeometry();
        EquationSystems* eq_sys = fe_data_manager->getEquationSystems();
        System& X_sys = eq_sys->get_system(fe_data_manager->COORDINATES_SYSTEM_NAME);
        FEDataManager::SystemDofMapCache& X_dof_map_cache =
            *fe_data_manager->getDofMapCache(fe_data_manager->COORDINATES_SYSTEM_NAME);
        NumericVector<double>* X_vec = X_sys.solution.get();

        IBTK::Point x_min, x_max;
        VectorValue<double> n;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const std::vector<Elem*>& active_patch_elems = active_patch_elem_map[patch->getPatchNumber()];
            if (active_patch_elems.size() == 0) continue;

            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
            Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(area_idx);
            Pointer<SideData<NDIM, double>> side_data = patch->getPatchData(side_idx);
            Pointer<NodeData<NDIM, double>> phi_data = patch->getPatchData(phi_idx);

            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const x_low = pgeom->getXLower();
            const double* const dx = pgeom->getDx();
            const hier::Index<NDIM>& idx_low = patch->getBox().lower();

            boost::multi_array<double, 2> x_node;
            using PointAxisSide = std::pair<libMesh::Point, std::pair<int, int>>;
            const IntVector<NDIM>& lower = patch->getBox().lower();
            const IntVector<NDIM>& upper = patch->getBox().upper();
            auto comp = [upper, lower](const hier::Index<NDIM>& a, const hier::Index<NDIM>& b) -> bool {
                int num_x = upper(0) - lower(0) + 1;
                int a_global = a(0) - lower(0) + num_x * (a(1) - lower(1) + 1);
                int b_global = b(0) - lower(0) + num_x * (b(1) - lower(1) + 1);
                return a_global < b_global;
            };
            std::map<hier::Index<NDIM>, std::vector<Elem*>, decltype(comp)> index_elem_map(comp);
            std::map<hier::Index<NDIM>, std::vector<PointAxisSide>, decltype(comp)> index_intersect_map(comp);
            for (const auto& elem : active_patch_elems)
            {
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                IBTK::get_values_for_interpolation(x_node, *X_vec, X_dof_indices);

                const unsigned int n_node = elem->n_nodes();
                std::vector<libMesh::Point> X_node_cache(n_node), x_node_cache(n_node);
                x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
                x_max = IBTK::Point::Constant(-std::numeric_limits<double>::max());
                for (unsigned int k = 0; k < n_node; ++k)
                {
                    X_node_cache[k] = elem->point(k);
                    libMesh::Point& x = x_node_cache[k];
                    for (unsigned int d = 0; d < NDIM; ++d) x(d) = x_node[k][d];

                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        x_min[d] = std::min(x_min[d], x(d));
                        x_max[d] = std::max(x_max[d], x(d));
                    }
                    elem->point(k) = x;
                }

                // Check if element is inside grid cell
                std::vector<hier::Index<NDIM>> elem_idx_nodes(n_node);
                for (unsigned int k = 0; k < n_node; ++k)
                {
                    const Node& node = elem->node_ref(k);
                    const hier::Index<NDIM>& idx = IndexUtilities::getCellIndex(&node(0), grid_geom, level->getRatio());
                    elem_idx_nodes[k] = idx;
                }
                if (std::adjacent_find(elem_idx_nodes.begin(), elem_idx_nodes.end(), std::not_equal_to<>()) ==
                    elem_idx_nodes.end())
                    TBOX_ERROR("Found an element completely contained within a grid cell.\n"
                               << "We are not currently equipped to handle these situations.\n");

                // Form bounding box of element
                Box<NDIM> box(IndexUtilities::getCellIndex(&x_min[0], grid_geom, level->getRatio()),
                              IndexUtilities::getCellIndex(&x_max[0], grid_geom, level->getRatio()));
                box.grow(1);
                box = box * (extended_box ? vol_data->getGhostBox() : patch->getBox());

                // We have the bounding box of the element. Loop over coordinate directions and look for intersections
                // with the background grid.
                for (BoxIterator<NDIM> b(box); b; b++)
                {
                    const hier::Index<NDIM>& i_c = b();
                    std::vector<libMesh::Point> intersection_points;
                    std::vector<std::pair<libMesh::Point, std::pair<int, int>>> axis_upper_lower_pairs;
                    bool added_elem = false;
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
                                r(d) = x_low[d] + dx[d] * (static_cast<double>(i_c(d) - idx_low(d)) +
                                                           (d == axis ? (upper_lower == 1 ? 1.0 : 0.0) : 0.5));
                            libMesh::Point p;

                            if (findIntersection(p, elem, r, q))
                            {
                                // We have an intersection
                                intersection_points.push_back(p);
                                // Add point to list of volume indices
                                if (!added_elem)
                                {
                                    index_elem_map[i_c].push_back(elem);
                                    added_elem = true;
                                }
                                // Make sure we haven't found this intersection before
                                const PointAxisSide pt_ax_si = std::make_pair(p, std::make_pair(axis, upper_lower));
                                auto pt_compare = [pt_ax_si](PointAxisSide a) -> bool {
                                    if (pt_ax_si.first == a.first)
                                        return true;
                                    else
                                        return false;
                                };
                                if (std::find_if(index_intersect_map[i_c].begin(),
                                                 index_intersect_map[i_c].end(),
                                                 pt_compare) == index_intersect_map[i_c].end())
                                    index_intersect_map[i_c].push_back(pt_ax_si);
                            }
                        }
                    }
                }
            }

            // We have all the intersections for this patch. Loop through and determine length fractions, surface areas,
            // and cell volumes
            for (const auto& index_elem_vec_pair : index_elem_map)
            {
                const hier::Index<NDIM>& idx = index_elem_vec_pair.first;
                const std::vector<Elem*>& elem_vec = index_elem_vec_pair.second;
                const std::vector<PointAxisSide>& pt_ax_si_vec = index_intersect_map[idx];
                if (pt_ax_si_vec.size() == 1)
                {
                    // There's only one intersection in this cell. This probably means the intersection is exactly on a
                    // node.
                    continue;
                }
                // Determine the "interior point" on the element, if there is one (There should be zero or one)
                libMesh::Point int_pt;
                bool has_int_pt = false;
                for (const auto& elem : elem_vec)
                {
                    for (unsigned int node_num = 0; node_num < elem->n_nodes(); ++node_num)
                    {
                        const libMesh::Point& pt = elem->point(node_num);
                        const hier::Index<NDIM>& pt_idx =
                            IndexUtilities::getCellIndex(&pt(0), grid_geom, level->getRatio());
                        if (pt_idx == idx)
                        {
                            int_pt = pt;
                            has_int_pt = true;
                        }
                    }
                }
                // Build new elem, this will be useful for computing distances/volumes
                std::unique_ptr<Elem> new_elem = Elem::build(EDGE2);
                std::array<std::unique_ptr<Node>, 2> new_nodes_for_elem;
                new_elem->set_id(0);
                for (int i = 0; i < 2; ++i)
                {
                    new_nodes_for_elem[i] = libmesh_make_unique<Node>(pt_ax_si_vec[i].first);
                    new_elem->set_node(i) = new_nodes_for_elem[i].get();
                }
                // Determine area contributions
                if (area_idx != IBTK::invalid_index)
                {
                    // If there's only one element, then the two intersections correspond to the single element.
                    double area = 0.0;
                    if (elem_vec.size() == 1)
                    {
                        area = (new_elem->point(0) - new_elem->point(1)).norm();
                    }
                    else
                    {
                        // Take care of edge case when element has point exactly on a side
                        if (!has_int_pt)
                        {
                            area += (pt_ax_si_vec[0].first - pt_ax_si_vec[1].first).norm();
                        }
                        else
                        {
                            for (unsigned int i = 0; i < pt_ax_si_vec.size(); ++i)
                            {
                                std::vector<libMesh::Point> pts = { pt_ax_si_vec[i].first, int_pt };
                                area += (pts[1] - pts[0]).norm();
                            }
                        }
                    }
                    (*area_data)(idx) = area;
                }
                // Determine sign of nodes
                std::array<std::array<double, 2>, 2> node_dist;
                // Determine normal for parent elements
                std::vector<Vector3d> elem_normals;
                Vector3d e3 = Vector3d::UnitZ() * (d_use_inside ? 1.0 : -1.0);
                for (const auto& elem : elem_vec)
                {
                    Vector3d v, w;
                    v << elem->point(0)(0), elem->point(0)(1), 0.0;
                    w << elem->point(1)(0), elem->point(1)(1), 0.0;
                    Vector3d n = (w - v).cross(e3);
                    n.normalize();
                    elem_normals.push_back(n);
                }
                // Determine distance to nodes
                for (int x = 0; x < 2; ++x)
                {
                    for (int y = 0; y < 2; ++y)
                    {
                        Vector3d P = Vector3d::Zero();
                        for (int d = 0; d < NDIM; ++d) P(d) = static_cast<double>(idx(d) + (d == 0 ? x : y));
                        // Project P onto element
                        Vector3d avg_proj, avg_unit_normal;
                        avg_proj.setZero();
                        avg_unit_normal.setZero();
                        double min_dist = std::numeric_limits<double>::max();
                        int num_min = 0;
                        for (unsigned int i = 0; i < elem_normals.size(); ++i)
                        {
                            const Elem* const elem = elem_vec[i];
                            const Vector3d& n = elem_normals[i];
                            Vector3d v, w;
                            v << (elem->point(0)(0) - x_low[0]) / dx[0], (elem->point(0)(1) - x_low[1]) / dx[1], 0.0;
                            w << (elem->point(1)(0) - x_low[0]) / dx[0], (elem->point(1)(1) - x_low[1]) / dx[1], 0.0;
                            const double t = std::max(0.0, std::min(1.0, (P - v).dot(w - v) / (v - w).squaredNorm()));
                            const Vector3d proj = v + t * (w - v);
                            const double dist = (proj - P).norm();
                            if (dist < min_dist)
                            {
                                min_dist = dist;
                                avg_proj = proj;
                                avg_unit_normal = n;
                                num_min = 1;
                            }
                            else if (MathUtilities<double>::equalEps(dist, min_dist))
                            {
                                avg_proj += proj;
                                avg_unit_normal += n;
                                ++num_min;
                            }
                        }
                        avg_proj /= static_cast<double>(num_min);
                        avg_unit_normal /= static_cast<double>(num_min);
                        avg_unit_normal.normalize();

                        const double dist = (P - avg_proj).norm();
                        Vector3d phys_vec;
                        for (unsigned int d = 0; d < NDIM; ++d) phys_vec(d) = dx[d] * (P - avg_proj)[d];
                        double dist_phys = phys_vec.norm();
                        VectorNd X;
                        for (int d = 0; d < NDIM; ++d) X[d] = x_low[d] + dx[d] * (P(d) - idx_low(d));
                        //                        double sgn = (avg_unit_normal.dot(P - avg_proj) <= 0.0 ? -1.0 : 1.0);
                        double ls_val = std::max(dist_up(X, d_y_up), dist_low(X, d_y_low));
                        double sgn = ls_val <= 0.0 ? -1.0 : +1.0;
                        node_dist[x][y] = dist * sgn;
                        NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y));
                        (*phi_data)(n_idx) = dist_phys * sgn;
                    }
                }
                // We have signs of nodes, now we can compute length fractions
                if (side_idx != IBTK::invalid_index)
                {
                    for (const auto& pt_ax_si : pt_ax_si_vec)
                    {
                        const libMesh::Point& pt = pt_ax_si.first;
                        const int& axis = pt_ax_si.second.first;
                        const int& upper_lower = pt_ax_si.second.second;
                        VectorNd pt_vec_1, pt_vec_2;
                        pt_vec_1 << (pt(0) - x_low[0]) / dx[0], (pt(1) - x_low[1]) / dx[1];
                        // Determine which side fraction to compute
                        int x = (axis == 0 ? upper_lower : (node_dist[0][upper_lower] < 0.0 ? 0 : 1));
                        int y = (axis == 1 ? upper_lower : (node_dist[upper_lower][0] < 0.0 ? 0 : 1));
                        pt_vec_2[0] = static_cast<double>(idx(0) - idx_low(0) + x);
                        pt_vec_2[1] = static_cast<double>(idx(1) - idx_low(1) + y);
                        SideIndex<NDIM> si(idx, axis, upper_lower);
                        (*side_data)(si) = (pt_vec_2 - pt_vec_1).norm();
                    }
                }
                // And we can compute volumes
                if (vol_idx != IBTK::invalid_index)
                {
                    // If there is a element with a node interior to a cell, calculate that volume first
                    double vol = 0.0;
#if (0)
                    // This doesn't need to be used if level set is generated from original mesh.
                    // The level set sees the "average" distance, which should be good enough for volume.
                    if (has_int_pt)
                    {
                        std::array<libMesh::Point, 3> tri_pts = { pt_ax_si_vec[0].first,
                                                                  pt_ax_si_vec[1].first,
                                                                  int_pt };
                        std::array<Vector3d, 3> tri_vecs;
                        for (int ii = 0; ii < 3; ++ii)
                        {
                            for (int d = 0; d < NDIM; ++d)
                            {
                                tri_vecs[ii](d) = tri_pts[ii](d);
                            }
                            tri_vecs[ii](2) = 0.0;
                        }
                        vol += 0.5 * ((tri_vecs[1] - tri_vecs[0]).cross(tri_vecs[1] - tri_vecs[2])).norm();
                    }
#endif
                    int num_neg = 0;
                    for (int x = 0; x < 2; ++x)
                    {
                        for (int y = 0; y < 2; ++y)
                        {
                            if (node_dist[x][y] <= 0.0) ++num_neg;
                        }
                    }
                    // ElemCutter does not treat intersections near nodes correctly!!
                    // We can compute volumes by sudividing our domain into triangles.
                    if (num_neg == 1)
                    {
                        // If there is only a single "negative node", then we have a triangle.
                        VectorNd pt0, pt1, pt2;

                        pt0 << ((pt_ax_si_vec[0].first)(0) - x_low[0]) / dx[0],
                            ((pt_ax_si_vec[0].first)(1) - x_low[1]) / dx[1];
                        pt1 << ((pt_ax_si_vec[1].first)(0) - x_low[0]) / dx[0],
                            ((pt_ax_si_vec[1].first)(1) - x_low[1]) / dx[1];
                        for (int x = 0; x < 2; ++x)
                        {
                            for (int y = 0; y < 2; ++y)
                            {
                                if (node_dist[x][y] <= 0.0)
                                {
                                    pt2 << static_cast<double>(idx(0) + x), static_cast<double>(idx(1) + y);
                                }
                            }
                        }
                        vol += (pt2 - pt0).norm() * (pt2 - pt1).norm() * 0.5;
                    }
                    else if (num_neg == 2)
                    {
                        // If there are two "negative nodes", then draw a line from one intersection to a negative node
                        // to get two triangles.
                        std::array<std::array<VectorNd, 3>, 2> simplices;
                        VectorNd pt0, pt1, pt2, pt3;
                        pt0 << ((pt_ax_si_vec[0].first)(0) - x_low[0]) / dx[0],
                            ((pt_ax_si_vec[0].first)(1) - x_low[1]) / dx[1];
                        pt1 << ((pt_ax_si_vec[1].first)(0) - x_low[0]) / dx[0],
                            ((pt_ax_si_vec[1].first)(1) - x_low[1]) / dx[1];
                        pt3.setZero();
                        for (int x = 0; x < 2; ++x)
                        {
                            for (int y = 0; y < 2; ++y)
                            {
                                if (node_dist[x][y] <= 0.0)
                                {
                                    // pt3 is on the line with pt0 that is perpendicular to a coordinate axis.
                                    VectorNd pt;
                                    pt << static_cast<double>(idx(0) + x), static_cast<double>(idx(1) + y);
                                    if (std::abs((pt - pt0).dot(VectorNd::UnitX())) < sqrt(DBL_EPSILON) ||
                                        std::abs((pt - pt0).dot(VectorNd::UnitY())) < sqrt(DBL_EPSILON))
                                        pt3 = pt;
                                    else
                                        pt2 = pt;
                                }
                            }
                        }
                        simplices[0] = { pt0, pt3, pt2 };
                        simplices[1] = { pt1, pt0, pt2 };
                        for (const auto& simplex : simplices)
                        {
                            const VectorNd &pt1 = simplex[0], pt2 = simplex[1], pt3 = simplex[2];
                            double a = (pt1 - pt2).norm(), b = (pt2 - pt3).norm(), c = (pt1 - pt3).norm();
                            double p = 0.5 * (a + b + c);
                            vol += std::sqrt(p * (p - a) * (p - b) * (p - c));
                        }
                    }
                    else if (num_neg == 3)
                    {
                        // If there are three "negative nodes", then the "positive" nodes form a triangle, and the
                        // volume is 1 - vol(triangle).
                        VectorNd pt0, pt1, pt2;
                        pt0 << ((pt_ax_si_vec[0].first)(0) - x_low[0]) / dx[0],
                            ((pt_ax_si_vec[0].first)(1) - x_low[1]) / dx[1];
                        pt1 << ((pt_ax_si_vec[1].first)(0) - x_low[0]) / dx[0],
                            ((pt_ax_si_vec[1].first)(1) - x_low[1]) / dx[1];
                        for (int x = 0; x < 2; ++x)
                        {
                            for (int y = 0; y < 2; ++y)
                            {
                                if (node_dist[x][y] > 0.0)
                                {
                                    pt2 << static_cast<double>(idx(0) + x), static_cast<double>(idx(1) + y);
                                }
                            }
                        }
                        vol += 1.0 - (pt2 - pt0).norm() * (pt2 - pt1).norm() * 0.5;
                    }
                    else if (num_neg == 4)
                    {
                        vol = 1.0;
                    }
                    (*vol_data)(idx) = vol;
                }
            }
        }
    }
    // Loop through the rest of the (ghosted) grid cells and calculate the distance based on the projection to the
    // initial location (NOTE: These values don't have much significance) We will have to take into account cut cells
    // that are outside the physical domain. We need "approximate" values for this (extend structure outside of domain).
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
        Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(area_idx);
        Pointer<SideData<NDIM, double>> side_data = patch->getPatchData(side_idx);
        Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(phi_idx);

        // Also loop over boundary boxes
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const xlow = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        const hier::Index<NDIM>& idx_low = patch->getBox().lower();
        std::vector<Box<NDIM>> fill_boxes;
        for (int d = 1; d <= NDIM; ++d)
        {
            const tbox::Array<BoundaryBox<NDIM>>& bdry_boxes = pgeom->getCodimensionBoundaries(d);
            for (int i = 0; i < bdry_boxes.size(); ++i)
            {
                const BoundaryBox<NDIM>& bdry_box = bdry_boxes[i];
                const int location_index = bdry_box.getLocationIndex();
                const int axis = location_index % 2;
                const int upper_lower = location_index / 2;
                if (pgeom->getTouchesRegularBoundary(axis, upper_lower))
                    fill_boxes.push_back(
                        pgeom->getBoundaryFillBox(bdry_box, patch->getBox(), vol_data->getGhostCellWidth()));
            }
        }
        for (const auto& box : fill_boxes)
        {
            for (CellIterator<NDIM> ci(box); ci; ci++)
            {
                const CellIndex<NDIM> idx = ci();
                double& vol = (*vol_data)(idx);
                if (vol == s_large_val)
                {
                    // We need to change this value.
                    VectorNd X_loc;
                    for (int x = 0; x < 2; ++x)
                    {
                        X_loc[0] = xlow[0] + dx[0] * static_cast<double>(idx(0) - idx_low(0) + x);
                        for (int y = 0; y < 2; ++y)
                        {
                            X_loc[1] = xlow[1] + dx[1] * static_cast<double>(idx(1) - idx_low(1) + y);
                            NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y));
                            double ls_val = std::max(dist_up(X_loc, d_y_up), dist_low(X_loc, d_y_low));
                            (*ls_data)(n_idx) = ls_val;
                        }
                    }
                    findVolume(xlow, dx, idx_low, ls_data, idx, vol);
                    vol /= (dx[0] * dx[1]);
                    for (int f = 0; f < 2; ++f)
                    {
#if (NDIM == 2)
                        double L = length_fraction(1.0,
                                                   (*ls_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(f, 0))),
                                                   (*ls_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(f, 1))));
#endif
                        (*side_data)(SideIndex<NDIM>(idx, 0, f)) = L;
                    }
                    for (int f = 0; f < 2; ++f)
                    {
#if (NDIM == 2)
                        double L = length_fraction(1.0,
                                                   (*ls_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(0, f))),
                                                   (*ls_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(1, f))));
#endif
                        (*side_data)(SideIndex<NDIM>(idx, 1, f)) = L;
                    }
                }
            }
        }
        // Calculate volumes based on the interior or exterior
        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            VectorNd X_loc;
            for (int d = 0; d < NDIM; ++d) X_loc[d] = xlow[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d) + 0.5);
            const double ls_val = std::max(dist_up(X_loc, d_y_up), dist_low(X_loc, d_y_low));
            double& vol = (*vol_data)(idx);
            if ((vol == s_large_val || vol == 1.0 || vol == 0.0 || vol != vol) && ls_val > 0.0)
            {
                vol = 0.0;
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (int upperlower = 0; upperlower < 2; ++upperlower)
                    {
                        SideIndex<NDIM> sidx(idx, axis, upperlower);
                        (*side_data)(sidx) = 0.0;
                    }
                }
                VectorNd X_loc;
                for (int x = 0; x < 2; ++x)
                {
                    X_loc[0] = xlow[0] + dx[0] * static_cast<double>(idx(0) - idx_low(0) + x);
                    for (int y = 0; y < 2; ++y)
                    {
                        X_loc[1] = xlow[1] + dx[1] * static_cast<double>(idx(1) - idx_low(1) + y);
                        NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y));
                        double ls_val = std::max(dist_up(X_loc, d_y_up), dist_low(X_loc, d_y_low));
                        (*ls_data)(n_idx) = ls_val;
                    }
                }
            }
            else if ((vol == s_large_val || vol == 1.0 || vol == 0.0 || vol != vol) && ls_val < 0.0)
            {
                vol = 1.0;
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (int upperlower = 0; upperlower < 2; ++upperlower)
                    {
                        SideIndex<NDIM> sidx(idx, axis, upperlower);
                        (*side_data)(sidx) = 1.0;
                    }
                }
                VectorNd X_loc;
                for (int x = 0; x < 2; ++x)
                {
                    X_loc[0] = xlow[0] + dx[0] * static_cast<double>(idx(0) - idx_low(0) + x);
                    for (int y = 0; y < 2; ++y)
                    {
                        X_loc[1] = xlow[1] + dx[1] * static_cast<double>(idx(1) - idx_low(1) + y);
                        NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y));
                        double ls_val = std::max(dist_up(X_loc, d_y_up), dist_low(X_loc, d_y_low));
                        (*ls_data)(n_idx) = ls_val;
                    }
                }
            }
        }
    }
    plog << "Checking volume data for NaNs\n";
    DebuggingUtilities::checkCellDataForNaNs(vol_idx, d_hierarchy, true);
    plog << "\nChecking area data for NaNs\n";
    DebuggingUtilities::checkCellDataForNaNs(area_idx, d_hierarchy, true);
    plog << "\nChecking side data for NaNs\n";
    DebuggingUtilities::checkSideDataForNaNs(side_idx, d_hierarchy, true);
    plog << "\nChecking LS data for NaNs\n";
    DebuggingUtilities::checkNodeDataForNaNs(phi_idx, d_hierarchy, true);
    LS_TIMER_STOP(t_updateVolumeAreaSideLS);
}

bool
LSPipeFlow::findIntersection(libMesh::Point& p, Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q)
{
    LS_TIMER_START(t_findIntersection);
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
    LS_TIMER_STOP(t_findIntersection);
    return found_intersection;
}

void
LSPipeFlow::findVolume(const double* const xlow,
                       const double* const dx,
                       const hier::Index<NDIM>& patch_lower,
                       Pointer<NodeData<NDIM, double>> phi_data,
                       const CellIndex<NDIM>& idx,
                       double& volume)
{
    // Create the initial simplices.
    std::vector<Simplex> simplices;
    // Create a vector of pairs of points and phi values
    VectorNd X;
    double phi;
    int num_p = 0, num_n = 0;
#if (NDIM == 2)
    boost::multi_array<std::pair<VectorNd, double>, NDIM> indices(boost::extents[2][2]);
#endif
#if (NDIM == 3)
    boost::multi_array<std::pair<VectorNd, double>, NDIM> indices(boost::extents[2][2][2]);
#endif
#if (NDIM == 2)
    for (int x = 0; x <= 1; ++x)
    {
        X(0) = xlow[0] + dx[0] * (idx(0) - patch_lower(0) + x);
        for (int y = 0; y <= 1; ++y)
        {
            X(1) = xlow[1] + dx[1] * (idx(1) - patch_lower(1) + y);
            NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y));
            phi = (*phi_data)(n_idx);
            if (std::abs(phi) < s_eps) phi = phi < 0.0 ? -s_eps : s_eps;
            indices[x][y] = std::make_pair(X, phi);
            if (phi > 0)
            {
                // Found a positive phi
                num_p++;
            }
            else
            {
                // Found a negative phi
                num_n++;
            }
        }
    }
#endif
#if (NDIM == 3)
    for (int x = 0; x <= 1; ++x)
    {
        X(0) = xlow[0] + dx[0] * (idx(0) - patch_lower(0) + x);
        for (int y = 0; y <= 1; ++y)
        {
            X(1) = xlow[1] + dx[1] * (idx(1) - patch_lower(1) + y);
            for (int z = 0; z <= 1; ++z)
            {
                X(2) = xlow[2] + dx[2] * (idx(2) - patch_lower(2) + z);
                NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y, z));
                phi = (*phi_data)(n_idx);
                indices[x][y][z] = std::make_pair(X, phi);
                if (phi > 0)
                {
                    num_p++;
                }
                else
                {
                    num_n++;
                }
            }
        }
    }
#endif
#if (NDIM == 2)
    // Divide grid cell in half to form two simplices.
    simplices.push_back({ indices[0][0], indices[1][0], indices[1][1] });
    simplices.push_back({ indices[0][0], indices[0][1], indices[1][1] });
#endif
#if (NDIM == 3)
    // Divide grid cell to form simplices.
    simplices.push_back({ indices[0][0][0], indices[1][0][0], indices[0][1][0], indices[0][0][1] });
    simplices.push_back({ indices[1][1][0], indices[1][0][0], indices[0][1][0], indices[1][1][1] });
    simplices.push_back({ indices[1][0][1], indices[1][0][0], indices[1][1][1], indices[0][0][1] });
    simplices.push_back({ indices[0][1][1], indices[1][1][1], indices[0][1][0], indices[0][0][1] });
    simplices.push_back({ indices[1][1][1], indices[1][0][0], indices[0][1][0], indices[0][0][1] });
#endif
    if (num_n == NDIM * NDIM)
    {
        // Grid cell is completely contained within physical boundary.
        volume = dx[0] * dx[1];
    }
    else if (num_p == NDIM * NDIM)
    {
        // Grid cell is completely outside of physical boundary.
        volume = 0.0;
    }
    else
    {
        volume = findVolume(simplices);
    }
}

double
LSPipeFlow::findVolume(const std::vector<Simplex>& simplices)
{
    // Loop over simplices
    std::vector<std::array<VectorNd, NDIM + 1>> final_simplices;
    for (const auto& simplex : simplices)
    {
        std::vector<int> n_phi, p_phi;
        for (size_t k = 0; k < simplex.size(); ++k)
        {
            const std::pair<VectorNd, double>& pt_pair = simplex[k];
            double phi = pt_pair.second;
            if (phi < 0)
            {
                n_phi.push_back(k);
            }
            else
            {
                p_phi.push_back(k);
            }
        }
        // Determine new simplices
#if (NDIM == 2)
        VectorNd pt0, pt1, pt2;
        double phi0, phi1, phi2;
        if (n_phi.size() == 1)
        {
            pt0 = simplex[n_phi[0]].first;
            pt1 = simplex[p_phi[0]].first;
            pt2 = simplex[p_phi[1]].first;
            phi0 = simplex[n_phi[0]].second;
            phi1 = simplex[p_phi[0]].second;
            phi2 = simplex[p_phi[1]].second;
            // Simplex is between P0, P01, P02
            VectorNd P01 = midpoint_value(pt0, phi0, pt1, phi1);
            VectorNd P02 = midpoint_value(pt0, phi0, pt2, phi2);
            final_simplices.push_back({ pt0, P01, P02 });
        }
        else if (n_phi.size() == 2)
        {
            pt0 = simplex[n_phi[0]].first;
            pt1 = simplex[n_phi[1]].first;
            pt2 = simplex[p_phi[0]].first;
            phi0 = simplex[n_phi[0]].second;
            phi1 = simplex[n_phi[1]].second;
            phi2 = simplex[p_phi[0]].second;
            // Simplex is between P0, P1, P02
            VectorNd P02 = midpoint_value(pt0, phi0, pt2, phi2);
            final_simplices.push_back({ pt0, pt1, P02 });
            // and P1, P12, P02
            VectorNd P12 = midpoint_value(pt1, phi1, pt2, phi2);
            final_simplices.push_back({ pt1, P12, P02 });
        }
        else if (n_phi.size() == 3)
        {
            pt0 = simplex[n_phi[0]].first;
            pt1 = simplex[n_phi[1]].first;
            pt2 = simplex[n_phi[2]].first;
            final_simplices.push_back({ pt0, pt1, pt2 });
        }
        else if (n_phi.size() == 0)
        {
            continue;
        }
        else
        {
            TBOX_ERROR("This statement should not be reached!");
        }
#endif
#if (NDIM == 3)
        VectorNd pt0, pt1, pt2, pt3;
        double phi0, phi1, phi2, phi3;
        if (n_phi.size() == 1)
        {
            pt0 = simplex[n_phi[0]].first;
            pt1 = simplex[p_phi[0]].first;
            pt2 = simplex[p_phi[1]].first;
            pt3 = simplex[p_phi[2]].first;
            phi0 = simplex[n_phi[0]].second;
            phi1 = simplex[p_phi[0]].second;
            phi2 = simplex[p_phi[1]].second;
            phi3 = simplex[p_phi[2]].second;
            // Simplex is between P0, P01, P02, P03
            VectorNd P01 = midpoint_value(pt0, phi0, pt1, phi1);
            VectorNd P02 = midpoint_value(pt0, phi0, pt2, phi2);
            VectorNd P03 = midpoint_value(pt0, phi0, pt3, phi3);
            final_simplices.push_back({ pt0, P01, P02, P03 });
        }
        else if (n_phi.size() == 2)
        {
            pt0 = simplex[n_phi[0]].first;
            pt1 = simplex[n_phi[1]].first;
            pt2 = simplex[p_phi[0]].first;
            pt3 = simplex[p_phi[1]].first;
            phi0 = simplex[n_phi[0]].second;
            phi1 = simplex[n_phi[1]].second;
            phi2 = simplex[p_phi[0]].second;
            phi3 = simplex[p_phi[1]].second;
            // Simplices are between P0, P1, P02, P13
            VectorNd P02 = midpoint_value(pt0, phi0, pt2, phi2);
            VectorNd P13 = midpoint_value(pt1, phi1, pt3, phi3);
            final_simplices.push_back({ pt0, pt1, P02, P13 });
            // and P12, P1, P02, P13
            VectorNd P12 = midpoint_value(pt1, phi1, pt2, phi2);
            final_simplices.push_back({ P12, pt1, P02, P13 });
            // and P0, P03, P02, P13
            VectorNd P03 = midpoint_value(pt0, phi0, pt3, phi3);
            final_simplices.push_back({ pt0, P03, P02, P13 });
        }
        else if (n_phi.size() == 3)
        {
            pt0 = simplex[n_phi[0]].first;
            pt1 = simplex[n_phi[1]].first;
            pt2 = simplex[n_phi[2]].first;
            pt3 = simplex[p_phi[0]].first;
            phi0 = simplex[n_phi[0]].second;
            phi1 = simplex[n_phi[1]].second;
            phi2 = simplex[n_phi[2]].second;
            phi3 = simplex[p_phi[0]].second;
            // Simplex is between P0, P1, P2, P13
            VectorNd P13 = midpoint_value(pt1, phi1, pt3, phi3);
            final_simplices.push_back({ pt0, pt1, pt2, P13 });
            // and P0, P03, P2, P13
            VectorNd P03 = midpoint_value(pt0, phi0, pt3, phi3);
            final_simplices.push_back({ pt0, P03, pt2, P13 });
            // and P23, P03, P2, P13
            VectorNd P23 = midpoint_value(pt2, phi2, pt3, phi3);
            final_simplices.push_back({ P23, P03, pt2, P13 });
        }
        else if (n_phi.size() == 4)
        {
            pt0 = simplex[n_phi[0]].first;
            pt1 = simplex[n_phi[1]].first;
            pt2 = simplex[n_phi[2]].first;
            pt3 = simplex[n_phi[3]].first;
            final_simplices.push_back({ pt0, pt1, pt2, pt3 });
        }
        else if (n_phi.size() == 0)
        {
            continue;
        }
        else
        {
            TBOX_ERROR("This statement should not be reached!");
        }
#endif
    }
    // Loop over simplices and compute volume
    double volume = 0.0;
    for (const auto& simplex : final_simplices)
    {
#if (NDIM == 2)
        VectorNd pt1 = simplex[0], pt2 = simplex[1], pt3 = simplex[2];
        double a = (pt1 - pt2).norm(), b = (pt2 - pt3).norm(), c = (pt1 - pt3).norm();
        double p = 0.5 * (a + b + c);
        volume += std::sqrt(p * (p - a) * (p - b) * (p - c));
#endif
#if (NDIM == 3)
        // Volume is given by 1/NDIM! * determinant of matrix
        Eigen::MatrixXd A(NDIM, NDIM);
        for (int d = 0; d < NDIM; ++d)
        {
            A.col(d) = simplex[d + 1] - simplex[0];
        }
        volume += 1.0 / 6.0 * std::abs(A.determinant());
#endif
    }
    return volume;
}
} // namespace LS
