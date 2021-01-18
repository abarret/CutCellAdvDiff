#include "LS/MLSReconstructCache.h"

#include "ibamr/app_namespaces.h"

namespace LS
{
MLSReconstructCache::MLSReconstructCache(int ls_idx,
                                         int vol_idx,
                                         Pointer<PatchHierarchy<NDIM>> hierarchy,
                                         bool use_centroids)
    : ReconstructCache(ls_idx, vol_idx, hierarchy, use_centroids)
{
    // intentionally blank
}

void
MLSReconstructCache::cacheData()
{
    // Can't cache data with MLS.
    d_update_weights = false;
    return;
}

double
MLSReconstructCache::reconstructOnIndex(VectorNd x_loc,
                                        const hier::Index<NDIM>& idx,
                                        const CellData<NDIM, double>& Q_data,
                                        Pointer<Patch<NDIM>> patch)
{
    Box<NDIM> box(idx, idx);
    box.grow(d_stencil_size);
    Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
    Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);
#ifndef NDEBUG
    TBOX_ASSERT(ls_data->getGhostBox().contains(box));
    TBOX_ASSERT(Q_data.getGhostBox().contains(box));
    TBOX_ASSERT(vol_data->getGhostBox().contains(box));
#endif
    int size = NDIM + 1;
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const xlow = pgeom->getXLower();
    const hier::Index<NDIM>& idx_low = patch->getBox().lower();
    for (int d = 0; d < NDIM; ++d) x_loc[d] = idx_low(d) + (x_loc[d] - xlow[d]) / dx[d];

    std::vector<double> Q_vals;
    std::vector<VectorNd> X_vals;

    for (CellIterator<NDIM> ci(box); ci; ci++)
    {
        const CellIndex<NDIM>& idx_c = ci();
        if ((*vol_data)(idx_c) > 0.0)
        {
            // Use this point to calculate least squares reconstruction.
            // Find cell center
            VectorNd x_cent_c;
            if (d_use_centroids)
            {
                x_cent_c = find_cell_centroid(idx_c, *ls_data);
            }
            else
            {
                for (int d = 0; d < NDIM; ++d) x_cent_c[d] = static_cast<double>(idx_c[d]) + 0.5;
            }
            Q_vals.push_back(Q_data(idx_c));
            X_vals.push_back(x_cent_c);
        }
    }
    const int m = Q_vals.size();
    VectorXd U(VectorXd::Zero(m));
    MatrixXd A(MatrixXd::Zero(m, size)), Lambda(MatrixXd::Zero(m, m));
    for (size_t i = 0; i < X_vals.size(); ++i)
    {
        const VectorNd X = X_vals[i] - x_loc;
        A(i, 2) = X[1];
        A(i, 1) = X[0];
        A(i, 0) = 1.0;
        Lambda(i, i) = std::sqrt(exp(static_cast<double>((X_vals[i] - x_loc).squaredNorm())));
        U(i) = Q_vals[i];
    }
    VectorXd x = (Lambda * A).fullPivHouseholderQr().solve(Lambda * U);
    return x(0);
}
} // namespace LS
