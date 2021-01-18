#ifndef included_LS_utility_functions
#define included_LS_utility_functions

#include "LS/LSFindCellVolume.h"
#include "LS/SetLSValue.h"

#include "Variable.h"
#include "tbox/MathUtilities.h"

#include "libmesh/elem.h"

#include "Eigen/Dense"

namespace LS
{
#define LS_TIMER_START(timer) timer->start();

#define LS_TIMER_STOP(timer) timer->stop();

static double s_eps = 1.0e-12;
inline double
length_fraction(const double dx, const double phi_l, const double phi_u)
{
    double L = 0.0;
    if (phi_l < 0.0 && phi_u > 0.0)
    {
        L = phi_l / (phi_l - phi_u);
    }
    else if (phi_l > 0.0 && phi_u < 0.0)
    {
        L = phi_u / (phi_u - phi_l);
    }
    else if (phi_l < 0.0 && phi_u < 0.0)
    {
        L = 1.0;
    }
    return L * dx;
}

inline double
area_fraction(const double reg_area, const double phi_ll, const double phi_lu, const double phi_uu, const double phi_ul)
{
    // Find list of vertices
    std::vector<IBTK::Vector2d> vertices;
    // Start at bottom left
    if (phi_ll < 0.0) vertices.push_back({ 0.0, 0.0 });
    // Go clockwise towards top
    if (phi_ll * phi_lu < 0.0) vertices.push_back({ 0.0, -phi_ll / (phi_lu - phi_ll) });
    if (phi_lu < 0.0) vertices.push_back({ 0.0, 1.0 });
    if (phi_lu * phi_uu < 0.0) vertices.push_back({ -phi_lu / (phi_uu - phi_lu), 1.0 });
    if (phi_uu < 0.0) vertices.push_back({ 1.0, 1.0 });
    if (phi_uu * phi_ul < 0.0) vertices.push_back({ 1.0, 1.0 - phi_uu / (phi_ul - phi_uu) });
    if (phi_ul < 0.0) vertices.push_back({ 1.0, 0.0 });
    if (phi_ul * phi_ll < 0.0) vertices.push_back({ 1.0 - phi_ul / (phi_ll - phi_ul), 0.0 });

    // We have vertices, now use shoelace formula to find area
    double A = 0.0;
    for (unsigned int i = 0; i < vertices.size(); ++i)
    {
        const IBTK::Vector2d& vertex = vertices[i];
        const IBTK::Vector2d& vertex_n = vertices[(i + 1) % vertices.size()];
        A += vertex(0) * vertex_n(1) - vertex_n(0) * vertex(1);
    }
    return 0.5 * std::abs(A) * reg_area;
}

inline IBTK::VectorNd
midpoint_value(const IBTK::VectorNd& pt0, const double& phi0, const IBTK::VectorNd& pt1, const double& phi1)
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
using Simplex = std::array<std::pair<VectorNd, double>, NDIM + 1>;
inline IBTK::VectorNd
find_cell_centroid_slow(const SAMRAI::pdat::CellIndex<NDIM>& idx, const SAMRAI::pdat::NodeData<NDIM, double>& ls_data)
{
    IBTK::VectorNd center;
    center.setZero();
    // Compute vertices
    VectorNd X;
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
inline VectorNd
find_cell_centroid(const SAMRAI::pdat::CellIndex<NDIM>& idx, const SAMRAI::pdat::NodeData<NDIM, double>& ls_data)
{
    IBTK::VectorNd centroid;
    std::vector<VectorNd> vertices;
    // Find all vertices
    for (int normal = 0; normal < NDIM; ++normal)
    {
        for (int side1 = 0; side1 < 2; ++side1)
        {
            for (int side2 = 0; side2 < 2; ++side2)
            {
                SAMRAI::hier::IntVector<NDIM> dir(0);
                dir(normal) = 1;
                IntVector<NDIM> low(0, 0, 0);
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

static inline double
rbf(double r)
{
    return r;
}

/*!
 * \brief Routine for converting strings to enums.
 */
template <typename T>
inline T
string_to_enum(const std::string& /*val*/)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return -1;
}

/*!
 * \brief Routine for converting enums to strings.
 */
template <typename T>
inline std::string enum_to_string(T /*val*/)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return "UNKNOWN";
}

enum class LeastSquaresOrder
{
    CONSTANT,
    LINEAR,
    QUADRATIC,
    CUBIC,
    UNKNOWN_ORDER = -1
};

template <>
inline LeastSquaresOrder
string_to_enum<LeastSquaresOrder>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "CONSTANT") == 0) return LeastSquaresOrder::CONSTANT;
    if (strcasecmp(val.c_str(), "LINEAR") == 0) return LeastSquaresOrder::LINEAR;
    if (strcasecmp(val.c_str(), "QUADRATIC") == 0) return LeastSquaresOrder::QUADRATIC;
    if (strcasecmp(val.c_str(), "CUBIC") == 0) return LeastSquaresOrder::CUBIC;
    return LeastSquaresOrder::UNKNOWN_ORDER;
}

template <>
inline std::string
enum_to_string<LeastSquaresOrder>(LeastSquaresOrder val)
{
    if (val == LeastSquaresOrder::CONSTANT) return "CONSTANT";
    if (val == LeastSquaresOrder::LINEAR) return "LINEAR";
    if (val == LeastSquaresOrder::QUADRATIC) return "QUADRATIC";
    if (val == LeastSquaresOrder::CUBIC) return "CUBIC";
    return "UNKNOWN_ORDER";
}

enum class AdvectionTimeIntegrationMethod
{
    FORWARD_EULER,
    MIDPOINT_RULE,
    UNKNOWN_METHOD
};

template <>
inline AdvectionTimeIntegrationMethod
string_to_enum<AdvectionTimeIntegrationMethod>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "FORWARD_EULER") == 0) return AdvectionTimeIntegrationMethod::FORWARD_EULER;
    if (strcasecmp(val.c_str(), "MIDPOINT_RULE") == 0) return AdvectionTimeIntegrationMethod::MIDPOINT_RULE;
    return AdvectionTimeIntegrationMethod::UNKNOWN_METHOD;
}

template <>
inline std::string
enum_to_string<AdvectionTimeIntegrationMethod>(AdvectionTimeIntegrationMethod val)
{
    if (val == AdvectionTimeIntegrationMethod::FORWARD_EULER) return "FORWARD_EULER";
    if (val == AdvectionTimeIntegrationMethod::MIDPOINT_RULE) return "MIDPOINT_RULE";
    return "UNKNOWN_METHOD";
}

enum class DiffusionTimeIntegrationMethod
{
    BACKWARD_EULER,
    TRAPEZOIDAL_RULE,
    UNKNOWN_METHOD
};

template <>
inline DiffusionTimeIntegrationMethod
string_to_enum<DiffusionTimeIntegrationMethod>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "BACKWARD_EULER") == 0) return DiffusionTimeIntegrationMethod::BACKWARD_EULER;
    if (strcasecmp(val.c_str(), "TRAPEZOIDAL_RULE") == 0) return DiffusionTimeIntegrationMethod::TRAPEZOIDAL_RULE;
    return DiffusionTimeIntegrationMethod::UNKNOWN_METHOD;
}

template <>
inline std::string
enum_to_string<DiffusionTimeIntegrationMethod>(DiffusionTimeIntegrationMethod val)
{
    if (val == DiffusionTimeIntegrationMethod::BACKWARD_EULER) return "BACKWARD_EULER";
    if (val == DiffusionTimeIntegrationMethod::TRAPEZOIDAL_RULE) return "TRAPEZOIDAL_RULE";
    return "UNKNOWN_METHOD";
}

enum class RBFPolyOrder
{
    LINEAR,
    QUADRATIC,
    UNKNOWN_ORDER
};

template <>
inline RBFPolyOrder
string_to_enum<RBFPolyOrder>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "LINEAR") == 0) return RBFPolyOrder::LINEAR;
    if (strcasecmp(val.c_str(), "QUADRATIC") == 0) return RBFPolyOrder::QUADRATIC;
    return RBFPolyOrder::UNKNOWN_ORDER;
}

template <>
inline std::string
enum_to_string<RBFPolyOrder>(RBFPolyOrder val)
{
    if (val == RBFPolyOrder::LINEAR) return "LINEAR";
    if (val == RBFPolyOrder::QUADRATIC) return "QUADRATIC";
    return "UNKNOWN_ORDER";
}

using ReactionFcn =
    std::function<double(double, const std::vector<double>&, const std::vector<double>&, double, void*)>;
using ReactionFcnCtx = std::pair<ReactionFcn, void*>;
using BdryConds = std::tuple<ReactionFcn, ReactionFcn, void*>;

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

struct PatchIndexPair
{
public:
    PatchIndexPair(const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch,
                   const SAMRAI::pdat::CellIndex<NDIM>& idx)
        : d_idx(idx)
    {
        d_patch_num = patch->getPatchNumber();
        const SAMRAI::hier::Box<NDIM>& box = patch->getBox();
        const SAMRAI::hier::Index<NDIM>& idx_low = box.lower();
        const SAMRAI::hier::Index<NDIM>& idx_up = box.upper();
        int num_x = idx_up(0) - idx_low(0) + 1;
        d_global_idx = idx(0) - idx_low(0) + num_x * (idx(1) - idx_low(1) + 1);
#if (NDIM == 3)
        int num_y = idx_up(1) - idx_low(1) + 1;
        d_global_idx += num_x * num_y * (idx(2) - idx_low(2));
#endif
    }

    bool operator<(const PatchIndexPair& b) const
    {
        bool less_than_b = false;
        if (d_patch_num < b.d_patch_num)
        {
            // We're on a smaller patch
            less_than_b = true;
        }
        else if (d_patch_num == b.d_patch_num && d_global_idx < b.d_global_idx)
        {
            // Our global index is smaller than b but on the same patch
            less_than_b = true;
        }
        return less_than_b;
    }

    SAMRAI::pdat::CellIndex<NDIM> d_idx;
    int d_patch_num = -1;
    int d_global_idx = -1;
};

struct CutCellElems
{
public:
    CutCellElems(libMesh::Elem* parent_elem, std::vector<libMesh::Point> intersections)
        : d_parent_elem(parent_elem), d_intersections(intersections)
    {
        d_elem = libMesh::Elem::build(d_parent_elem->type());
        for (size_t n = 0; n < d_intersections.size(); ++n)
        {
            d_nodes.push_back(libMesh::Node::build(d_intersections[n], n));
            d_elem->set_id(0);
            d_elem->set_node(n) = d_nodes[n].get();
        }
    }

    libMesh::Elem* d_parent_elem;
    std::unique_ptr<libMesh::Elem> d_elem;
    std::vector<std::unique_ptr<libMesh::Node>> d_nodes;
    std::vector<libMesh::Point> d_intersections;
};

} // namespace LS
#endif
