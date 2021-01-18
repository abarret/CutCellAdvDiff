#include "LS/IntegrateFunction.h"
#include "LS/utility_functions.h"

#include "tbox/ShutdownRegistry.h"

#include "boost/multi_array.hpp"

#include "ibamr/namespaces.h"

namespace LS
{
/////////////////// STATIC /////////////////////////////

IntegrateFunction* IntegrateFunction::s_integrate_function = nullptr;
unsigned char IntegrateFunction::s_shutdown_priority = 225;
#if (NDIM == 2)
const std::array<double, 9> IntegrateFunction::s_weights = { 0.9876542474e-1, 0.1391378575e-1, 0.1095430035,
                                                             0.6172839460e-1, 0.8696116674e-2, 0.6846438175e-1,
                                                             0.6172839460e-1, 0.8696116674e-2, 0.6846438175e-1 };
const std::array<IBTK::VectorNd, 9> IntegrateFunction::s_quad_pts = { IBTK::VectorNd(0.25, 0.5),
                                                                      IBTK::VectorNd(0.5635083269e-1, 0.8872983346),
                                                                      IBTK::VectorNd(0.4436491673, 0.1127016654),
                                                                      IBTK::VectorNd(0.4436491673, 0.5),
                                                                      IBTK::VectorNd(0.1, 0.8872983346),
                                                                      IBTK::VectorNd(0.7872983346, 0.1127016654),
                                                                      IBTK::VectorNd(0.5635083269e-1, 0.5),
                                                                      IBTK::VectorNd(0.1270166538e-1, 0.8872983346),
                                                                      IBTK::VectorNd(0.1, 0.1127016654) };
#endif
#if (NDIM == 3)
const std::array<double, 27> IntegrateFunction::s_weights = {
    3.068198819728420e-5, 3.068198819728420e-5, 4.909118111565470e-5, 2.177926162424280e-4, 2.177926162424280e10 - 4,
    3.484681859878840e-4, 2.415587821057510e-4, 2.415587821057510e-4, 3.864940513692010e-4, 9.662351284230000e-4,
    9.662351284230000e-4, 0.001545976205477,    0.006858710562414,    0.006858710562414,    0.010973936899863,
    0.007607153074595,    0.007607153074595,    0.012171444919352,    0.001901788268649,    0.001901788268649,
    0.003042861229838,    0.013499628508586,    0.013499628508586,    0.021599405613738,    0.014972747367084,
    0.014972747367084,    0.023956395787334
};

const std::array<IBTK::VectorNd, 27> IntegrateFunction::s_quad_pts = {
    VectorNd(0.001431498841332, 0.011270166537926, 0.100000000000000),
    VectorNd(0.001431498841332, 0.011270166537926, 0.100000000000000),
    VectorNd(0.006350832689629, 0.006350832689629, 0.100000000000000),
    VectorNd(0.006350832689629, 0.050000000000000, 0.056350832689629),
    VectorNd(0.050000000000000, 0.006350832689629, 0.056350832689629),
    VectorNd(0.028175416344815, 0.028175416344815, 0.056350832689629),
    VectorNd(0.011270166537926, 0.088729833462074, 0.012701665379258),
    VectorNd(0.088729833462074, 0.011270166537926, 0.012701665379258),
    VectorNd(0.050000000000000, 0.050000000000000, 0.012701665379258),
    VectorNd(0.006350832689629, 0.050000000000000, 0.443649167310371),
    VectorNd(0.050000000000000, 0.006350832689629, 0.443649167310371),
    VectorNd(0.028175416344815, 0.028175416344815, 0.443649167310371),
    VectorNd(0.028175416344815, 0.221824583655185, 0.250000000000000),
    VectorNd(0.221824583655185, 0.028175416344815, 0.250000000000000),
    VectorNd(0.125000000000000, 0.125000000000000, 0.250000000000000),
    VectorNd(0.050000000000000, 0.393649167310371, 0.056350832689629),
    VectorNd(0.393649167310371, 0.050000000000000, 0.056350832689629),
    VectorNd(0.221824583655185, 0.221824583655185, 0.056350832689629),
    VectorNd(0.011270166537926, 0.088729833462074, 0.787298334620741),
    VectorNd(0.088729833462074, 0.011270166537926, 0.787298334620741),
    VectorNd(0.050000000000000, 0.050000000000000, 0.787298334620741),
    VectorNd(0.050000000000000, 0.393649167310371, 0.443649167310371),
    VectorNd(0.393649167310371, 0.050000000000000, 0.443649167310371),
    VectorNd(0.221824583655185, 0.221824583655185, 0.443649167310371),
    VectorNd(0.088729833462074, 0.698568501158667, 0.100000000000000),
    VectorNd(0.698568501158667, 0.088729833462074, 0.100000000000000),
    VectorNd(0.393649167310371, 0.393649167310371, 0.100000000000000)
};
#endif
const double IntegrateFunction::s_eps = 1.0e-12;

IntegrateFunction*
IntegrateFunction::getIntegrator()
{
    if (!s_integrate_function)
    {
        s_integrate_function = new IntegrateFunction();
        ShutdownRegistry::registerShutdownRoutine(freeIntegrators, s_shutdown_priority);
    }
    return s_integrate_function;
}

void
IntegrateFunction::freeIntegrators()
{
    if (s_integrate_function) delete s_integrate_function;
    s_integrate_function = nullptr;
}

void
IntegrateFunction::integrateFcnOnPatchHierarchy(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                const int ls_idx,
                                                const int Q_idx,
                                                fcn_type fcn,
                                                double t)
{
    d_ls_idx = ls_idx;
    d_hierarchy = hierarchy;
    d_fcn = fcn;
    d_t = t;

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
            Q_data->fillAll(0.0);
            Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(ls_idx);

            const Box<NDIM>& box = patch->getBox();
            const hier::Index<NDIM>& idx_l = box.lower();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();

            for (CellIterator<NDIM> ci(box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();

                VectorNd XLowerCorner;
                for (int d = 0; d < NDIM; ++d)
                    XLowerCorner(d) = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_l(d)));
                (*Q_data)(idx) = integrateOverIndex(dx, XLowerCorner, *ls_data, idx);
            }
        }
    }
    return;
}

/////////////////// PRIVATE ////////////////////////////

double
IntegrateFunction::integrateOverIndex(const double* const dx,
                                      const VectorNd& XLow,
                                      const NodeData<NDIM, double>& ls_data,
                                      const CellIndex<NDIM>& idx)
{
    // Create initial simplices
    std::vector<Simplex> simplices;
    // Create a vector of pairs of points and phi values
    VectorNd X;
    double phi;
    int num_p = 0;
#if (NDIM == 2)
    boost::multi_array<std::pair<VectorNd, double>, NDIM> indices(boost::extents[2][2]);
    for (int x = 0; x <= 1; ++x)
    {
        X(0) = XLow(0) + dx[0] * x;
        for (int y = 0; y <= 1; ++y)
        {
            X(1) = XLow(1) + dx[1] * y;
            NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y));
            phi = ls_data(n_idx);
            if (std::abs(phi) < s_eps) phi = phi < 0.0 ? -s_eps : s_eps;
            indices[x][y] = std::make_pair(X, phi);
            if (phi > 0) ++num_p;
        }
    }
    simplices.push_back({ indices[0][0], indices[1][0], indices[1][1] });
    simplices.push_back({ indices[0][0], indices[0][1], indices[1][1] });
#endif
#if (NDIM == 3)
    boost::multi_array<std::pair<VectorNd, double>, NDIM> indices(boost::extents[2][2][2]);
    for (int x = 0; x <= 1; ++x)
    {
        X(0) = XLow(0) + dx[0] * x;
        for (int y = 0; y <= 1; ++y)
        {
            X(1) = XLow(1) + dx[1] * y;
            for (int z = 0; z <= 1; ++z)
            {
                X(2) = XLow(2) + dx[2] * z;
                NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y, z));
                phi = ls_data(n_idx);
                indices[x][y][z] = std::make_pair(X, phi);
                if (phi > 0) ++num_p;
            }
        }
    }
    simplices.push_back({ indices[0][0][0], indices[1][0][0], indices[0][1][0], indices[0][0][1] });
    simplices.push_back({ indices[1][1][0], indices[1][0][0], indices[0][1][0], indices[1][1][1] });
    simplices.push_back({ indices[1][0][1], indices[1][0][0], indices[1][1][1], indices[0][0][1] });
    simplices.push_back({ indices[0][1][1], indices[1][1][1], indices[0][1][0], indices[0][0][1] });
    simplices.push_back({ indices[1][1][1], indices[1][0][0], indices[0][1][0], indices[0][0][1] });
#endif
    double vol = 0.0;
    if (num_p == NDIM * NDIM)
        vol = 0.0;
    else
        vol = integrate(simplices);
    return vol;
}

double
IntegrateFunction::integrate(const std::vector<Simplex>& simplices)
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
            if (phi < 0.0)
                n_phi.push_back(k);
            else
                p_phi.push_back(k);
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
            VectorNd P02 = midpoint_value(pt0, phi0, pt2, phi2);
            final_simplices.push_back({ pt0, pt1, P02 });
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
    // Loop over simplices and compute integral
    double vol = 0.0;
    for (const auto& simplex : final_simplices)
    {
        vol += integrateOverSimplex(simplex, d_t);
    }
    return vol;
}

double
IntegrateFunction::integrateOverSimplex(const std::array<VectorNd, NDIM + 1>& X_pts, const double t)
{
#if (NDIM == 3)
    double integral = 0.0;
    for (const auto& X_pt : X_pts)
    {
        integral += d_fcn(X_pt, t);
    }
    MatrixXd mat = MatrixXd::Zero(NDIM, NDIM);
    for (size_t l = 1; l < X_pts.size(); ++l)
    {
        mat.block(0, l - 1, NDIM, 1) = X_pts[l] - X_pts[0];
    }
    integral *= abs(mat.determinant()) / 6.0;
    return integral;
#endif
#if (NDIM == 2)
    double integral = 0.0;
    double J = std::abs((X_pts[1](0) - X_pts[0](0)) * (X_pts[2](1) - X_pts[0](1)) -
                        (X_pts[2](0) - X_pts[0](0)) * (X_pts[1](1) - X_pts[0](1)));
#ifndef NDEBUG
    if (J <= 0.0)
    {
        pout << "pt 0: \n" << X_pts[0] << "\n";
        pout << "pt 1: \n" << X_pts[1] << "\n";
        pout << "pt 2: \n" << X_pts[2] << "\n";
        TBOX_ERROR("Found zero or negative J: " << J << "\n");
    }
#endif
    for (size_t i = 0; i < s_weights.size(); ++i)
    {
        // Convert reference coords to physical coords
        VectorNd X = referenceToPhysical(s_quad_pts[i], X_pts);
        // Evaluate function
        double val = d_fcn(X, t);
        if (val > 1.0e10) pout << "val is large.\n";
        integral += val * s_weights[i];
    }
    return integral * J;
#endif
}
} // namespace LS
