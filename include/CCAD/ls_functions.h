#ifndef included_CCAD_utility_functions
#define included_CCAD_utility_functions
#include "CCAD/LSFindCellVolume.h"
#include "CCAD/ls_utilities.h"

#include "ibtk/ibtk_utilities.h"

#include "NodeIndex.h"
#include "Variable.h"
#include "tbox/MathUtilities.h"

#include "libmesh/elem.h"
#include "libmesh/vector_value.h"

#include "boost/multi_array.hpp"

#include "Eigen/Dense"

namespace CCAD
{
#define CCAD_TIMER_START(timer) timer->start();

#define CCAD_TIMER_STOP(timer) timer->stop();

static double s_eps = 1.0e-12;

inline SAMRAI::pdat::NodeIndex<NDIM>
get_node_index_from_corner(const SAMRAI::hier::Index<NDIM>& idx, int corner)
{
    switch (corner)
    {
    case 0:
        return SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::pdat::NodeIndex<NDIM>::LowerLeft);
    case 1:
        return SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::pdat::NodeIndex<NDIM>::LowerRight);
    case 2:
        return SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::pdat::NodeIndex<NDIM>::UpperRight);
    case 3:
        return SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::pdat::NodeIndex<NDIM>::UpperLeft);
    default:
        TBOX_ERROR("Unknown corner number\n");
        return SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::pdat::NodeIndex<NDIM>::LowerLeft);
    }
}
double length_fraction(double dx, double phi_l, double phi_u);

double area_fraction(double reg_area, double phi_ll, double phi_lu, double phi_uu, double phi_ul);

inline IBTK::VectorNd
midpoint_value(const IBTK::VectorNd& pt0, double phi0, const IBTK::VectorNd& pt1, double phi1)
{
    return pt0 * phi1 / (phi1 - phi0) - pt1 * phi0 / (phi1 - phi0);
}

#if (NDIM == 2)
inline IBTK::VectorNd
find_cell_centroid(const SAMRAI::pdat::CellIndex<NDIM>& idx, const SAMRAI::pdat::NodeData<NDIM, double>& ls_data)
{
    IBTK::VectorNd center;
    center.setZero();
    std::vector<IBTK::VectorNd> X_pts;

    double phi_ll = ls_data(SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::hier::IntVector<NDIM>(0, 0)));
    if (std::abs(phi_ll) < s_eps) phi_ll = phi_ll < 0.0 ? -s_eps : s_eps;
    double phi_ul = ls_data(SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::hier::IntVector<NDIM>(1, 0)));
    if (std::abs(phi_ul) < s_eps) phi_ul = phi_ul < 0.0 ? -s_eps : s_eps;
    double phi_uu = ls_data(SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::hier::IntVector<NDIM>(1, 1)));
    if (std::abs(phi_uu) < s_eps) phi_uu = phi_uu < 0.0 ? -s_eps : s_eps;
    double phi_lu = ls_data(SAMRAI::pdat::NodeIndex<NDIM>(idx, SAMRAI::hier::IntVector<NDIM>(0, 1)));
    if (std::abs(phi_lu) < s_eps) phi_lu = phi_lu < 0.0 ? -s_eps : s_eps;
    if ((phi_ll < 0.0 && phi_ul < 0.0 && phi_uu < 0.0 && phi_lu < 0.0) ||
        (phi_ll > 0.0 && phi_ul > 0.0 && phi_uu > 0.0 && phi_lu > 0.0))
    {
        // Not a cut cell. Center is idx
        center(0) = idx(0) + 0.5;
        center(1) = idx(1) + 0.5;
    }
    else
    {
        // Loop over nodes and edges and find points
        if (phi_ll < 0.0) X_pts.push_back(IBTK::VectorNd(idx(0), idx(1)));
        if (phi_ll * phi_lu < 0.0) X_pts.push_back(IBTK::VectorNd(idx(0), idx(1) - phi_ll / (phi_lu - phi_ll)));
        if (phi_lu < 0.0) X_pts.push_back(IBTK::VectorNd(idx(0), idx(1) + 1.0));
        if (phi_lu * phi_uu < 0.0) X_pts.push_back(IBTK::VectorNd(idx(0) - phi_lu / (phi_uu - phi_lu), idx(1) + 1.0));
        if (phi_uu < 0.0) X_pts.push_back(IBTK::VectorNd(idx(0) + 1.0, idx(1) + 1.0));
        if (phi_uu * phi_ul < 0.0) X_pts.push_back(IBTK::VectorNd(idx(0) + 1.0, idx(1) - phi_ul / (phi_uu - phi_ul)));
        if (phi_ul < 0.0) X_pts.push_back(IBTK::VectorNd(idx(0) + 1.0, idx(1)));
        if (phi_ul * phi_ll < 0.0) X_pts.push_back(IBTK::VectorNd(idx(0) - phi_ll / (phi_ul - phi_ll), idx(1)));

        double signed_area = 0.0;
        for (size_t i = 0; i < X_pts.size(); ++i)
        {
            const IBTK::VectorNd& X = X_pts[i];
            const IBTK::VectorNd& X_n = X_pts[(i + 1) % X_pts.size()];
            center += (X + X_n) * (X(0) * X_n(1) - X_n(0) * X(1));
            signed_area += 0.5 * (X(0) * X_n(1) - X_n(0) * X(1));
        }
        center /= 6.0 * signed_area;
        if (std::abs(signed_area) < 1.0e-8)
        {
            // Degenerate polygon. Switch to average of vertices.
            center.setZero();
            for (const auto& X : X_pts)
            {
                center += X;
            }
            center /= static_cast<double>(X_pts.size());
        }
    }

    return center;
}
#endif
#if (NDIM == 3)
// Slow, accurate computation of cell centroid
using Simplex = std::array<std::pair<IBTK::VectorNd, double>, NDIM + 1>;
inline IBTK::VectorNd
find_cell_centroid_slow(const SAMRAI::pdat::CellIndex<NDIM>& idx, const SAMRAI::pdat::NodeData<NDIM, double>& ls_data)
{
    IBTK::VectorNd center;
    center.setZero();
    // Compute vertices
    IBTK::VectorNd X;
    double phi;
    int num_p = 0, num_n = 0;
    boost::multi_array<std::pair<IBTK::VectorNd, double>, NDIM> indices(boost::extents[2][2][2]);
    for (int x = 0; x <= 1; ++x)
    {
        X(0) = static_cast<double>(idx(0) + x);
        for (int y = 0; y <= 1; ++y)
        {
            X(1) = static_cast<double>(idx(1) + y);
            for (int z = 0; z <= 1; ++z)
            {
                X(2) = static_cast<double>(idx(2) + z);
                SAMRAI::pdat::NodeIndex<NDIM> n_idx(idx, SAMRAI::hier::IntVector<NDIM>(x, y, z));
                phi = ls_data(n_idx);
                if (std::abs(phi) < s_eps) phi = phi < 0.0 ? -s_eps : s_eps;
                indices[x][y][z] = std::make_pair(X, phi);
                if (phi > 0.0)
                    ++num_p;
                else
                    ++num_n;
            }
        }
    }
    if ((num_n == NDIM * NDIM) || (num_p == NDIM * NDIM))
    {
        // Grid cell is interior or exterior, cell centroid is idx + 0.5
        for (int d = 0; d < NDIM; ++d) center[d] = static_cast<double>(idx(d)) + 0.5;
        return center;
    }

    // We are on a cut cell.
    // Break apart vertices into simplices. Then integrate along the surface
    std::vector<Simplex> simplices;
    simplices.push_back({ indices[0][0][0], indices[1][0][0], indices[0][1][0], indices[0][0][1] });
    simplices.push_back({ indices[1][1][0], indices[1][0][0], indices[0][1][0], indices[1][1][1] });
    simplices.push_back({ indices[1][0][1], indices[1][0][0], indices[1][1][1], indices[0][0][1] });
    simplices.push_back({ indices[0][1][1], indices[1][1][1], indices[0][1][0], indices[0][0][1] });
    simplices.push_back({ indices[1][1][1], indices[1][0][0], indices[0][1][0], indices[0][0][1] });

    std::vector<std::array<IBTK::VectorNd, NDIM + 1>> final_simplices;
    for (const auto& simplex : simplices)
    {
        // Loop through simplices.
        // Determine boundary points on simplex to form polytope.
        std::vector<int> n_phi, p_phi;
        for (unsigned long k = 0; k < simplex.size(); ++k)
        {
            const std::pair<IBTK::VectorNd, double>& pt_pair = simplex[k];
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
        IBTK::VectorNd pt0, pt1, pt2, pt3;
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
            IBTK::VectorNd P01 = midpoint_value(pt0, phi0, pt1, phi1);
            IBTK::VectorNd P02 = midpoint_value(pt0, phi0, pt2, phi2);
            IBTK::VectorNd P03 = midpoint_value(pt0, phi0, pt3, phi3);
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
            IBTK::VectorNd P02 = midpoint_value(pt0, phi0, pt2, phi2);
            IBTK::VectorNd P13 = midpoint_value(pt1, phi1, pt3, phi3);
            final_simplices.push_back({ pt0, pt1, P02, P13 });
            // and P12, P1, P02, P13
            IBTK::VectorNd P12 = midpoint_value(pt1, phi1, pt2, phi2);
            final_simplices.push_back({ P12, pt1, P02, P13 });
            // and P0, P03, P02, P13
            IBTK::VectorNd P03 = midpoint_value(pt0, phi0, pt3, phi3);
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
            IBTK::VectorNd P13 = midpoint_value(pt1, phi1, pt3, phi3);
            final_simplices.push_back({ pt0, pt1, pt2, P13 });
            // and P0, P03, P2, P13
            IBTK::VectorNd P03 = midpoint_value(pt0, phi0, pt3, phi3);
            final_simplices.push_back({ pt0, P03, pt2, P13 });
            // and P23, P03, P2, P13
            IBTK::VectorNd P23 = midpoint_value(pt2, phi2, pt3, phi3);
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
    }
    // Loop through final_simplices and calculate centroid of each tetrahedron
    for (const auto& simplex : final_simplices)
    {
        // Centroid of tetrahedron is average of vertices
        IBTK::VectorNd tet_cent;
        tet_cent.setZero();
        for (const auto& vertex : simplex) tet_cent += vertex;
        tet_cent /= static_cast<double>(simplex.size());
        center += tet_cent;
    }
    for (int d = 0; d < NDIM; ++d) center(d) /= static_cast<double>(final_simplices.size());
    return center;
}
// Fast, but possibly wrong computation of cell centroid
inline IBTK::VectorNd
find_cell_centroid(const SAMRAI::pdat::CellIndex<NDIM>& idx, const SAMRAI::pdat::NodeData<NDIM, double>& ls_data)
{
    IBTK::VectorNd centroid;
    std::vector<IBTK::VectorNd> vertices;
    // Find all vertices
    for (int normal = 0; normal < NDIM; ++normal)
    {
        for (int side1 = 0; side1 < 2; ++side1)
        {
            for (int side2 = 0; side2 < 2; ++side2)
            {
                SAMRAI::hier::IntVector<NDIM> dir(0);
                dir(normal) = 1;
                SAMRAI::hier::IntVector<NDIM> low(0, 0, 0);
                if (normal == 0)
                {
                    low(1) = side1;
                    low(2) = side2;
                }
                else if (normal == 1)
                {
                    low(0) = side1;
                    low(2) = side2;
                }
                else
                {
                    low(0) = side1;
                    low(1) = side2;
                }
                const SAMRAI::pdat::NodeIndex<NDIM> ni(idx, low);
                const double& phi_l = ls_data(ni);
                const double& phi_u = ls_data(ni + dir);

                IBTK::VectorNd x_l = IBTK::VectorNd::Zero(), x_u = IBTK::VectorNd::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    x_l(d) = static_cast<double>(ni(d));
                    x_u(d) = static_cast<double>(ni(d) + dir(d));
                }
                if (phi_l * phi_u < 0.0) vertices.push_back(midpoint_value(x_l, phi_l, x_u, phi_u));
            }
        }
    }
    // Loop through vertices of box
    for (int x = 0; x < 2; ++x)
    {
        for (int y = 0; y < 2; ++y)
        {
            for (int z = 0; z < 2; ++z)
            {
                SAMRAI::hier::IntVector<NDIM> node(x, y, z);
                const SAMRAI::pdat::NodeIndex<NDIM> ni(idx, node);
                if (ls_data(ni) < 0.0)
                    vertices.push_back(
                        { static_cast<double>(ni(0)), static_cast<double>(ni(1)), static_cast<double>(ni(2)) });
            }
        }
    }
    if (vertices.size() == 0)
    {
        for (int d = 0; d < NDIM; ++d) centroid[d] = static_cast<double>(idx(d)) + 0.5;
        return centroid;
    }

    centroid.setZero();
    for (const auto& vertex : vertices)
    {
        centroid += vertex;
    }
    for (int d = 0; d < NDIM; ++d) centroid(d) /= static_cast<double>(vertices.size());
    return centroid;
}
#endif

inline double
node_to_cell(const SAMRAI::pdat::CellIndex<NDIM>& idx, SAMRAI::pdat::NodeData<NDIM, double>& ls_data)
{
#if (NDIM == 2)
    SAMRAI::pdat::NodeIndex<NDIM> idx_ll(idx, SAMRAI::hier::IntVector<NDIM>(0, 0));
    SAMRAI::pdat::NodeIndex<NDIM> idx_lu(idx, SAMRAI::hier::IntVector<NDIM>(0, 1));
    SAMRAI::pdat::NodeIndex<NDIM> idx_ul(idx, SAMRAI::hier::IntVector<NDIM>(1, 0));
    SAMRAI::pdat::NodeIndex<NDIM> idx_uu(idx, SAMRAI::hier::IntVector<NDIM>(1, 1));
    double ls_ll = ls_data(idx_ll), ls_lu = ls_data(idx_lu);
    double ls_ul = ls_data(idx_ul), ls_uu = ls_data(idx_uu);
    return 0.25 * (ls_ll + ls_lu + ls_ul + ls_uu);
#endif
#if (NDIM == 3)
    double val = 0.0;
    for (int x = 0; x < 2; ++x)
    {
        for (int y = 0; y < 2; ++y)
        {
            for (int z = 0; z < 2; ++z)
            {
                SAMRAI::pdat::NodeIndex<NDIM> n_idx(idx, SAMRAI::hier::IntVector<NDIM>(x, y, z));
                val += ls_data(n_idx);
            }
        }
    }
    return val / 6.0;
#endif
}

inline void
copy_face_to_side(const int u_s_idx,
                  const int u_f_idx,
                  SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> s_data = patch->getPatchData(u_s_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double>> f_data = patch->getPatchData(u_f_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SAMRAI::pdat::SideIndex<NDIM>& s_idx = si();
                    SAMRAI::pdat::FaceIndex<NDIM> f_idx(s_idx.toCell(0), axis, 1);
                    (*s_data)(s_idx) = (*f_data)(f_idx);
                }
            }
        }
    }
}

static inline bool
findIntersection(libMesh::Point& p, libMesh::Elem* elem, const libMesh::Point& r, const libMesh::VectorValue<double>& q)
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

static inline std::string
get_libmesh_restart_file_name(const std::string& restart_dump_dirname,
                              const std::string& base_filename,
                              unsigned int time_step_number,
                              unsigned int part,
                              const std::string& extension)
{
    std::ostringstream file_name_prefix;
    file_name_prefix << restart_dump_dirname << "/" << base_filename << "_part_" << part << "." << std::setw(6)
                     << std::setfill('0') << std::right << time_step_number << "." << extension;
    return file_name_prefix.str();
}
} // namespace CCAD
#endif
