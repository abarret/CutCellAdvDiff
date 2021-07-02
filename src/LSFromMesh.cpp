#include "CCAD/LSFromMesh.h"
#include "CCAD/ls_functions.h"

#include "ibtk/DebuggingUtilities.h"
#include "ibtk/IBTK_MPI.h"

#include "RefineAlgorithm.h"

#if (0)
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Quotient.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_triangulation_2.h>
#endif

#include <algorithm>
#include <queue>

namespace
{
static Timer* t_updateVolumeAreaSideLS = nullptr;
} // namespace

namespace CCAD
{
LSFromMesh::LSFromMesh(std::string object_name,
                       Pointer<PatchHierarchy<NDIM>> hierarchy,
                       const Pointer<CutCellMeshMapping>& cut_cell_mesh_mapping,
                       bool use_inside /* = true*/)
    : LSFindCellVolume(std::move(object_name), hierarchy),
      d_use_inside(use_inside),
      d_cut_cell_mesh_mapping(cut_cell_mesh_mapping),
      d_sgn_var(new NodeVariable<NDIM, double>(d_object_name + "::SGN_VAR"))
{
    commonConstructor();
    return;
} // Constructor

void
LSFromMesh::updateVolumeAreaSideLS(int vol_idx,
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
    CCAD_TIMER_START(t_updateVolumeAreaSideLS);
    TBOX_ASSERT(phi_idx != IBTK::invalid_index);
    TBOX_ASSERT(phi_var);
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM>>& patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> phi_data = patch->getPatchData(phi_idx);
            phi_data->fillAll(static_cast<double>(ln + 2));
        }
    }
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, 0, finest_ln);
    hier_cc_data_ops.setToScalar(vol_idx, 0.0, false);
    hier_cc_data_ops.setToScalar(area_idx, 0.0, false);
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(d_hierarchy, 0, finest_ln);
    hier_sc_data_ops.setToScalar(side_idx, 0.0, false);

    d_cut_cell_mesh_mapping->initializeObjectState(d_hierarchy);
    d_cut_cell_mesh_mapping->generateCutCellMappings();
    const std::vector<std::map<IndexList, std::vector<CutCellElems>>>& idx_cut_cell_map_vec =
        d_cut_cell_mesh_mapping->getIdxCutCellElemsMap(finest_ln);

    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(finest_ln);
    unsigned int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const std::map<IndexList, std::vector<CutCellElems>>& idx_cut_cell_map = idx_cut_cell_map_vec[local_patch_num];
        if (idx_cut_cell_map.size() == 0) continue;

        Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
        Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(area_idx);
        Pointer<SideData<NDIM, double>> side_data = patch->getPatchData(side_idx);
        Pointer<NodeData<NDIM, double>> phi_data = patch->getPatchData(phi_idx);

        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const x_low = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        const hier::Index<NDIM>& idx_low = patch->getBox().lower();

        // Loop through all cut cells in map
        for (const auto& idx_elem_vec_pair : idx_cut_cell_map)
        {
            const CellIndex<NDIM>& idx = idx_elem_vec_pair.first.d_idx;
            const std::vector<CutCellElems>& cut_cell_elem_vec = idx_elem_vec_pair.second;

            // Determine area contributions
            if (area_data)
            {
                double area = 0.0;
                for (const auto& cut_cell_elem : cut_cell_elem_vec)
                {
                    area += (cut_cell_elem.d_elem->point(1) - cut_cell_elem.d_elem->point(0)).norm();
                }
                (*area_data)(idx) = area;
            }

            // Let's find distances from background cell nodes to structure
            // Determine normal for elements
            // Warning, we need the normal to be consistent between parent and child elements.
            std::vector<IBTK::Vector3d> elem_normals;
            for (const auto& cut_cell_elem : cut_cell_elem_vec)
            {
                // Note we use the parent element to calculate normals to preserve directions
                Vector3d v, w;
                const std::array<libMesh::Point, 2>& parent_pts = cut_cell_elem.d_parent_cur_pts;
                const unsigned int part = cut_cell_elem.d_part;
                v << parent_pts[0](0), parent_pts[0](1), parent_pts[0](2);
                w << parent_pts[1](0), parent_pts[1](1), parent_pts[1](2);
                const unsigned int domain_id = cut_cell_elem.d_parent_elem->subdomain_id();
                Vector3d e3 = Vector3d::UnitZ();
                if (!d_use_inside) e3 *= -1.0;
                if (d_norm_reverse_domain_ids[part].find(domain_id) != d_norm_reverse_domain_ids[part].end() ||
                    d_norm_reverse_elem_ids[part].find(cut_cell_elem.d_parent_elem->id()) !=
                        d_norm_reverse_elem_ids[part].end())
                {
                    e3 *= -1.0;
                }
                Vector3d n = (w - v).cross(e3);
                elem_normals.push_back(n);
            }
            // Determine distances to nodes
            for (int x = 0; x < 2; ++x)
            {
                for (int y = 0; y < 2; ++y)
                {
                    Vector3d P = Vector3d::Zero();
                    for (int d = 0; d < NDIM; ++d) P(d) = static_cast<double>(idx(d) - idx_low(d)) + (d == 0 ? x : y);
                    // Project P onto element
                    Vector3d avg_proj, avg_unit_normal;
                    avg_proj.setZero();
                    avg_unit_normal.setZero();
                    double min_dist = std::numeric_limits<double>::max();
                    int num_min = 0;
                    // Loop through all elements and calculate the smallest distance
                    for (unsigned int i = 0; i < elem_normals.size(); ++i)
                    {
                        const std::unique_ptr<Elem>& elem = cut_cell_elem_vec[i].d_elem;
                        const Vector3d& n = elem_normals[i];
                        Vector3d v, w;
                        v << (elem->point(0)(0) - x_low[0]) / dx[0], (elem->point(0)(1) - x_low[1]) / dx[1], 0.0;
                        w << (elem->point(1)(0) - x_low[0]) / dx[0], (elem->point(1)(1) - x_low[1]) / dx[1], 0.0;
                        const double t = std::max(0.0, std::min(1.0, (P - v).dot(w - v) / (v - w).squaredNorm()));
                        const Vector3d proj = v + t * (w - v);
                        VectorNd x_proj;
                        for (int d = 0; d < NDIM; ++d) x_proj[d] = x_low[d] + dx[d] * (proj(d) - idx_low(d));
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

                    Vector3d phys_vec;
                    for (unsigned int d = 0; d < NDIM; ++d) phys_vec(d) = dx[d] * (P - avg_proj)[d];
                    double dist_phys = phys_vec.norm();
                    double sgn = (avg_unit_normal.dot(P - avg_proj) <= 0.0 ? -1.0 : 1.0);
                    NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y));
                    (*phi_data)(n_idx) =
                        dist_phys < std::abs((*phi_data)(n_idx)) ? (dist_phys * sgn) : (*phi_data)(n_idx);
                }
            }

#if (0)
            // TODO: Need to determine the correct way to account for regions with interior points.
            // Find cell volumes. We use a triangulation algorithm from CGAL
            // First find all vertices and put them in one vector
            using K = CGAL::Exact_predicates_inexact_constructions_kernel;
            using Triangulation = CGAL::Triangulation_2<K>;
            bool draw = false;
            if (idx(0) == 128 && idx(1) == 67) draw = true;
            Triangulation t;
            t.clear();
            std::vector<K::Point_2> points;
            // First vertices of cell nodes
            VectorNd x_pt;
            for (int d = 0; d < NDIM; ++d) x_pt[d] = x_low[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
            NodeIndex<NDIM> idx_ll(idx, NodeIndex<NDIM>::LowerLeft);
            NodeIndex<NDIM> idx_lr(idx, NodeIndex<NDIM>::LowerRight);
            NodeIndex<NDIM> idx_ur(idx, NodeIndex<NDIM>::UpperRight);
            NodeIndex<NDIM> idx_ul(idx, NodeIndex<NDIM>::UpperLeft);
            if ((*phi_data)(idx_ll) < 0.0) points.push_back(K::Point_2(x_pt[0], x_pt[1]));
            if ((*phi_data)(idx_lr) < 0.0) points.push_back(K::Point_2(x_pt[0] + dx[0], x_pt[1]));
            if ((*phi_data)(idx_ur) < 0.0) points.push_back(K::Point_2(x_pt[0] + dx[0], x_pt[1] + dx[1]));
            if ((*phi_data)(idx_ul) < 0.0) points.push_back(K::Point_2(x_pt[0], x_pt[1] + dx[1]));
            // Now loop through all intersections
            // Only count interior points once
            std::vector<libMesh::Point> int_pts;
            std::map<int, std::vector<libMesh::Point>> side_pt_map;
            for (const auto& cut_cell : cut_cell_elem_vec)
            {
                const std::vector<std::pair<libMesh::Point, int>>& pt_side_vec = cut_cell.d_intersection_side_vec;
                for (const auto& pt_side : pt_side_vec)
                {
                    const libMesh::Point& pt = pt_side.first;
                    if (pt_side.second == -1 && std::find(int_pts.begin(), int_pts.end(), pt) != int_pts.end())
                        continue; // Already found point.
                    else
                        int_pts.push_back(pt); // Haven't found point yet.
                    side_pt_map[pt_side.second].push_back(pt);
                    points.push_back(K::Point_2(pt(0), pt(1)));
                }
            }
            // If we have multiple intersections on a side, we need to check whether that side has INTERIOR cell nodes.
            // If there are INTERIOR cell nodes, our structure is concave, and we need to use a different algorithm
            bool concave_structure = false;
            for (const auto& side_pt_vec_pair : side_pt_map)
            {
                int side = side_pt_vec_pair.first;
                if (side == -1) continue;
                if (side_pt_vec_pair.second.size() > 1)
                {
                    NodeIndex<NDIM> idx_l, idx_u;
                    switch (side)
                    {
                    case 0:
                        idx_l = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::LowerLeft);
                        idx_u = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::UpperLeft);
                        break;
                    case 1:
                        idx_l = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::LowerRight);
                        idx_u = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::UpperRight);
                        break;
                    case 2:
                        idx_l = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::LowerLeft);
                        idx_u = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::LowerRight);
                        break;
                    case 3:
                        idx_l = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::UpperLeft);
                        idx_u = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::UpperRight);
                        break;
                    default:
                        TBOX_ERROR("Unknown side!\n");
                    }
                    if ((*phi_data)(idx_l) < 0.0 && (*phi_data)(idx_u) < 0.0) concave_structure = true;
                }
            }
            // We now have all the points in the vector. Use the triangulation algorithm.
            t.insert(points.begin(), points.end());
            if (draw) CGAL::draw(t);
            Triangulation::Finite_faces_iterator face_iter = t.finite_faces_begin();
            double vol = 0.0;
            for (; face_iter != t.finite_faces_end(); ++face_iter)
            {
                Triangulation::Triangle tri = t.triangle(face_iter);
                vol += tri.area();
            }
            (*vol_data)(idx) = vol / (dx[0] * dx[1]);
#else
            // Find volumes from level set.
            double vol = 0.0;
            findVolume(x_low, dx, idx_low, phi_data, idx, vol);
            (*vol_data)(idx) = vol / (dx[0] * dx[1]);
#endif

            // Determine side lengths
            if (side_idx != IBTK::invalid_index)
            {
                for (int f = 0; f < 2; ++f)
                {
#if (NDIM == 2)
                    double L = length_fraction(1.0,
                                               (*phi_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(f, 0))),
                                               (*phi_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(f, 1))));
#endif
#if (NDIM == 3)
                    double L = 0.0;
                    TBOX_ERROR("3D Not implemented yet.\n");
#endif
                    (*side_data)(SideIndex<NDIM>(idx, 0, f)) = L;
                }
                for (int f = 0; f < 2; ++f)
                {
#if (NDIM == 2)
                    double L = length_fraction(1.0,
                                               (*phi_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(0, f))),
                                               (*phi_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(1, f))));
#endif
#if (NDIM == 3)
                    double L = 0.0;
                    TBOX_ERROR("3D Not implemented yet.\n");
#endif
                    (*side_data)(SideIndex<NDIM>(idx, 1, f)) = L;
                }
#if (NDIM == 3)
                for (int f = 0; f < 2; ++f)
                {
                    double L = 0.0;
                    TBOX_ERROR("3D Not implemented yet.\n");
                    (*side_data)(SideIndex<NDIM>(idx, 2, f)) = L;
                }
#endif
            }
        }
    }

    // Fill in physical boundary cells
    if (d_bdry_fcn)
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(phi_idx);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
            Pointer<SideData<NDIM, double>> side_data = patch->getPatchData(side_idx);

            // Loop over boundary boxes
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
                            pgeom->getBoundaryFillBox(bdry_box, patch->getBox(), ls_data->getGhostCellWidth()));
                }
            }
            for (const auto& box : fill_boxes)
            {
                for (NodeIterator<NDIM> ni(box); ni; ni++)
                {
                    const NodeIndex<NDIM>& idx = ni();
                    if ((*ls_data)(idx) == static_cast<double>(finest_ln + 2))
                    {
                        // Change this value
                        VectorNd X_loc;
                        for (int d = 0; d < NDIM; ++d)
                            X_loc[d] = xlow[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
                        double ls_val = (*ls_data)(idx);
                        d_bdry_fcn(X_loc, ls_val);
                        (*ls_data)(idx) = ls_val;
                    }
                }
                // Now fill in Cell values
                for (CellIterator<NDIM> ci(box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    double& vol = (*vol_data)(idx);
                    if (vol > 1.0)
                    {
                        // We need to change this value.
                        findVolume(xlow, dx, idx_low, ls_data, idx, vol);
                        vol /= (dx[0] * dx[1]);
                        for (int f = 0; f < 2; ++f)
                        {
#if (NDIM == 2)
                            double L = length_fraction(1.0,
                                                       (*ls_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(f, 0))),
                                                       (*ls_data)(NodeIndex<NDIM>(idx, IntVector<NDIM>(f, 1))));
#endif
#if (NDIM == 3)
                            double L = 0.0;
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
#if (NDIM == 3)
                            double L = 0.0;
#endif
                            (*side_data)(SideIndex<NDIM>(idx, 1, f)) = L;
                        }
#if (NDIM == 3)
                        for (int f = 0; f < 2; ++f)
                        {
                            TBOX_ERROR("Three dimensions not supported yet.\n");
                            double L = 0.0;
                            (*side_data)(SideIndex<NDIM>(idx, 2, f)) = L;
                        }
#endif
                    }
                }
            }
        }
    }

    updateLSAwayFromInterface(phi_idx);
    // Finally, fill in volumes/areas of non cut cells
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
            Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(area_idx);
            Pointer<NodeData<NDIM, double>> sgn_data = patch->getPatchData(phi_idx);
            Pointer<SideData<NDIM, double>> side_data = patch->getPatchData(side_idx);

            const Box<NDIM>& box = vol_data->getGhostBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();

            for (CellIterator<NDIM> ci(box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                double phi_cc = node_to_cell(idx, *sgn_data);
                if (phi_cc > 0.0)
                {
                    if (vol_data && (*vol_data)(idx) == 0.0) (*vol_data)(idx) = 0.0;
                    if (area_data && (*area_data)(idx) == 0.0) (*area_data)(idx) = 0.0;
                    if (side_data)
                    {
                        for (int axis = 0; axis < NDIM; ++axis)
                        {
                            for (int upper_lower = 0; upper_lower < 2; ++upper_lower)
                            {
                                SideIndex<NDIM> si(idx, axis, upper_lower);
                                if ((*side_data)(si) == 0.0) (*side_data)(si) = 0.0;
                            }
                        }
                    }
                }
                else if (phi_cc < 0.0)
                {
                    if (vol_data && (*vol_data)(idx) == 0.0) (*vol_data)(idx) = 1.0;
                    if (area_data && (*area_data)(idx) == 0.0) (*area_data)(idx) = 0.0;
                    if (side_data)
                    {
                        for (int axis = 0; axis < NDIM; ++axis)
                        {
                            for (int upper_lower = 0; upper_lower < 2; ++upper_lower)
                            {
                                SideIndex<NDIM> si(idx, axis, upper_lower);
                                if ((*side_data)(si) == 0.0) (*side_data)(si) = 1.0;
                            }
                        }
                    }
                }
            }
        }
    }
    CCAD_TIMER_STOP(t_updateVolumeAreaSideLS);
}

void
LSFromMesh::commonConstructor()
{
    IBAMR_DO_ONCE(t_updateVolumeAreaSideLS =
                      TimerManager::getManager()->getTimer("CCAD::LSFromMesH::updateVolumeAreaSideLS()"););
    d_norm_reverse_domain_ids.resize(d_cut_cell_mesh_mapping->getNumParts());
    d_norm_reverse_elem_ids.resize(d_cut_cell_mesh_mapping->getNumParts());

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_sgn_idx = var_db->registerVariableAndContext(d_sgn_var, var_db->getContext(d_object_name + "::Context"), 1);
}

void
LSFromMesh::findVolume(const double* const xlow,
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
            if (std::abs(phi) < 1.0e-8) phi = phi < 0.0 ? -1.0e-8 : 1.0e-8;
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
LSFromMesh::findVolume(const std::vector<Simplex>& simplices)
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

void
LSFromMesh::updateLSAwayFromInterface(const int phi_idx)
{
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    HierarchyNodeDataOpsReal<NDIM, double> hier_nc_data_ops(d_hierarchy, 0, finest_ln);
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(finest_ln);
    // Now we need to update the sign of phi_data.
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_sgn_idx);
    }
    hier_nc_data_ops.copyData(d_sgn_idx, phi_idx, false);
    floodFillForLS(finest_ln, static_cast<double>(finest_ln + 2));
    // At this point, the finest level has been filled in. We now need to fill in coarser levels
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    Pointer<CoarsenOperator<NDIM>> coarsen_op = grid_geom->lookupCoarsenOperator(d_sgn_var, "CONSTANT_COARSEN");
    Pointer<CoarsenAlgorithm<NDIM>> coarsen_alg = new CoarsenAlgorithm<NDIM>();
    coarsen_alg->registerCoarsen(d_sgn_idx, d_sgn_idx, coarsen_op);
    std::vector<Pointer<CoarsenSchedule<NDIM>>> coarsen_scheds(finest_ln + 1);
    for (int ln = finest_ln; ln > 0; --ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM>> coarser_level = d_hierarchy->getPatchLevel(ln - 1);
        coarsen_scheds[ln] = coarsen_alg->createSchedule(coarser_level, level);
        coarsen_scheds[ln]->coarsenData();
        floodFillForLS(ln - 1, static_cast<double>(ln + 1));
    }

    hier_nc_data_ops.copyData(phi_idx, d_sgn_idx, false);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_sgn_idx);
    }
}

void
LSFromMesh::floodFillForLS(const int ln, const double eps)
{
    TBOX_ASSERT(eps > 0.0);
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(d_sgn_idx, d_sgn_idx, d_sgn_idx, nullptr);
    Pointer<RefineSchedule<NDIM>> ghost_fill_sched = ghost_fill_alg.createSchedule(level);
    // Do a flood fill algorithm
    std::vector<int> patch_filled_vec(level->getNumberOfPatches());
    unsigned int patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& box = patch->getBox();
        Box<NDIM> n_box(box);
        n_box.growUpper(1);
        Pointer<NodeData<NDIM, double>> phi_data = patch->getPatchData(d_sgn_idx);
        NodeData<NDIM, int> idx_touched(box, 1, phi_data->getGhostCellWidth());
        idx_touched.fillAll(0);
        std::queue<NodeIndex<NDIM>> idx_queue;
        bool found_pt = false;
        for (NodeIterator<NDIM> ni(box); ni; ni++)
        {
            const NodeIndex<NDIM>& idx = ni();
            if ((*phi_data)(idx) < 0.0)
            {
                idx_queue.push(idx);
                found_pt = true;
            }
        }
        patch_filled_vec[patch_num] = found_pt ? 1 : 0;
        // We have our starting point. Now, loop through queue
        while (idx_queue.size() > 0)
        {
            const NodeIndex<NDIM>& idx = idx_queue.front();
            // If this point is uninitialized, it is interior
            if (idx_touched(idx) == 0 && ((*phi_data)(idx) == eps || (*phi_data)(idx) < 0.0))
            {
                // Insert the point into touched list
                idx_touched(idx) = 1;
                if ((*phi_data)(idx) == eps) (*phi_data)(idx) = -eps;
                // Add neighboring points if they haven't been touched yet
                NodeIndex<NDIM> idx_s = idx + IntVector<NDIM>(0, -1);
                NodeIndex<NDIM> idx_n = idx + IntVector<NDIM>(0, 1);
                NodeIndex<NDIM> idx_e = idx + IntVector<NDIM>(1, 0);
                NodeIndex<NDIM> idx_w = idx + IntVector<NDIM>(-1, 0);
                if (n_box.contains(idx_s) && ((*phi_data)(idx_s) == eps || (*phi_data)(idx_s) < 0.0) &&
                    idx_touched(idx_s) == 0)
                    idx_queue.push(idx_s);
                if (n_box.contains(idx_n) && ((*phi_data)(idx_n) == eps || (*phi_data)(idx_n) < 0.0) &&
                    idx_touched(idx_n) == 0)
                    idx_queue.push(idx_n);
                if (n_box.contains(idx_e) && ((*phi_data)(idx_e) == eps || (*phi_data)(idx_e) < 0.0) &&
                    idx_touched(idx_e) == 0)
                    idx_queue.push(idx_e);
                if (n_box.contains(idx_w) && ((*phi_data)(idx_w) == eps || (*phi_data)(idx_w) < 0.0) &&
                    idx_touched(idx_w) == 0)
                    idx_queue.push(idx_w);
            }
            idx_queue.pop();
        }
    }

    // At this point if there's any box that hasn't been filled, then it's either entirely inside or outside.
    // We'll fill ghost cells, then check the ghost boxes. If there is any negative in the ghost box, then the entire
    // patch is inside
    int num_negative_found = 1;
    while (num_negative_found > 0)
    {
        num_negative_found = 0;
        ghost_fill_sched->fillData(0.0);
        patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
        {
            bool found_negative = false;
            if (patch_filled_vec[patch_num] == 1) continue;
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            // Loop through ghost cells
            Pointer<NodeData<NDIM, double>> phi_data = patch->getPatchData(d_sgn_idx);
            for (NodeIterator<NDIM> ni(phi_data->getGhostBox()); ni; ni++)
            {
                const NodeIndex<NDIM>& idx = ni();
                if (patch->getBox().contains(idx)) continue;
                if ((*phi_data)(idx) == -eps)
                {
                    found_negative = true;
                    break;
                }
            }
            if (found_negative)
            {
                phi_data->fillAll(-eps, phi_data->getGhostBox());
                num_negative_found++;
                patch_filled_vec[patch_num] = 1;
            }
        }
        num_negative_found = IBTK_MPI::sumReduction(num_negative_found);
    }
}
} // namespace CCAD
