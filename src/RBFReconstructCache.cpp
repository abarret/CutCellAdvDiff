#include "LS/RBFReconstructCache.h"

#include "ibamr/app_namespaces.h"

namespace LS
{
RBFReconstructCache::RBFReconstructCache(int ls_idx,
                                         int vol_idx,
                                         Pointer<PatchHierarchy<NDIM>> hierarchy,
                                         bool use_centroids)
    : ReconstructCache(ls_idx, vol_idx, hierarchy, use_centroids)
{
    // intentionally blank
}

void
RBFReconstructCache::cacheData()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    // Free any preallocated matrices
    for (std::map<PatchIndexPair, FullPivHouseholderQR<MatrixXd>>& qr_matrix_map : d_qr_matrix_vec)
        qr_matrix_map.clear();

    // allocate matrix data
    d_qr_matrix_vec.resize(finest_ln + 1);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::map<PatchIndexPair, FullPivHouseholderQR<MatrixXd>>& qr_map = d_qr_matrix_vec[ln];
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            const Box<NDIM>& box = patch->getBox();
            Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);

            for (CellIterator<NDIM> ci(box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                if ((*vol_data)(idx) < 1.0 && (*vol_data)(idx) > 0.0)
                {
                    // We are on a cut cell. We need to interpolate to cell center
                    VectorNd x_loc;
                    for (int d = 0; d < NDIM; ++d) x_loc(d) = static_cast<double>(idx(d)) + 0.5;
                    Box<NDIM> box(idx, idx);
                    box.grow(d_stencil_size);
#ifndef NDEBUG
                    TBOX_ASSERT(ls_data->getGhostBox().contains(box));
                    TBOX_ASSERT(vol_data->getGhostBox().contains(box));
#endif

                    const CellIndex<NDIM>& idx_low = patch->getBox().lower();
                    std::vector<VectorNd> X_vals;

                    for (CellIterator<NDIM> ci(box); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx_c = ci();
                        if ((*vol_data)(idx_c) > 0.0)
                        {
                            // Use this point to calculate least squares reconstruction.
                            VectorNd x_cent_c;
                            if (d_use_centroids)
                            {
                                x_cent_c = find_cell_centroid(idx_c, *ls_data);
                            }
                            else
                            {
                                for (int d = 0; d < NDIM; ++d) x_cent_c[d] = static_cast<double>(idx_c[d]) + 0.5;
                            }
                            X_vals.push_back(x_cent_c);
                        }
                    }
                    const int m = X_vals.size();
                    MatrixXd A(MatrixXd::Zero(m, m));
                    MatrixXd B(MatrixXd::Zero(m, NDIM + 1));
                    for (size_t i = 0; i < X_vals.size(); ++i)
                    {
                        for (size_t j = 0; j < X_vals.size(); ++j)
                        {
                            const VectorNd X = X_vals[i] - X_vals[j];
                            const double phi = rbf(X.norm());
                            A(i, j) = phi;
                        }
                        B(i, 0) = 1.0;
                        for (int d = 0; d < NDIM; ++d) B(i, d + 1) = X_vals[i](d);
                    }
                    PatchIndexPair p_idx = PatchIndexPair(patch, idx);

                    MatrixXd final_mat(MatrixXd::Zero(m + NDIM + 1, m + NDIM + 1));
                    final_mat.block(0, 0, m, m) = A;
                    final_mat.block(0, m, m, NDIM + 1) = B;
                    final_mat.block(m, 0, NDIM + 1, m) = B.transpose();

#ifndef NDEBUG
                    if (qr_map.find(p_idx) == qr_map.end())
                        qr_map[p_idx] = FullPivHouseholderQR<MatrixXd>(final_mat);
                    else
                        TBOX_WARNING("Already had a QR decomposition in place");
#else
                    qr_map[p_idx] = FullPivHouseholderQR<MatrixXd>(final_mat);
#endif
                }
            }
        }
    }
    d_update_weights = false;
}

double
RBFReconstructCache::reconstructOnIndex(VectorNd x_loc,
                                        const hier::Index<NDIM>& idx,
                                        const CellData<NDIM, double>& Q_data,
                                        Pointer<Patch<NDIM>> patch)
{
    if (d_update_weights) cacheData();
    const int ln = patch->getPatchLevelNumber();
    Box<NDIM> box(idx, idx);
    box.grow(d_stencil_size);
    Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
    Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);
#ifndef NDEBUG
    TBOX_ASSERT(ls_data->getGhostBox().contains(box));
    TBOX_ASSERT(Q_data.getGhostBox().contains(box));
    TBOX_ASSERT(vol_data->getGhostBox().contains(box));
#endif
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
    VectorXd U(VectorXd::Zero(m + NDIM + 1));
    for (size_t i = 0; i < Q_vals.size(); ++i) U(i) = Q_vals[i];

    PatchIndexPair pi_pair(patch, idx);
    VectorXd x1 = d_qr_matrix_vec[ln][PatchIndexPair(patch, idx)].solve(U);
    VectorXd rbf_coefs = x1.block(0, 0, m, 1);
    VectorXd poly_coefs = x1.block(m, 0, NDIM + 1, 1);
    VectorXd poly_vec = VectorXd::Ones(NDIM + 1);
    for (int d = 0; d < NDIM; ++d) poly_vec(d + 1) = x_loc(d);
    double val = 0.0;
    for (size_t i = 0; i < X_vals.size(); ++i)
    {
        val += rbf_coefs[i] * rbf((X_vals[i] - x_loc).norm());
    }
    val += poly_coefs.dot(poly_vec);
    return val;
}
} // namespace LS
