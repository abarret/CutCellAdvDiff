#ifndef included_ReconstructCache
#define included_ReconstructCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"

#include "LS/utility_functions.h"

#include "CellData.h"
#include "CellIndex.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"

#include <Eigen/Dense>

#include <map>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace LS
{
/*!
 * \brief Class ReconstructCache caches the data necessary to form RBF reconstructions of data.
 */
class ReconstructCache : public SAMRAI::tbox::DescribedClass
{
public:
    ReconstructCache() = default;

    ReconstructCache(int ls_idx,
                     int vol_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                     bool use_centroids = true);

    virtual ~ReconstructCache() = default;

    /*!
     * \brief Deleted copy constructor.
     */
    ReconstructCache(const ReconstructCache& from) = delete;

    virtual void clearCache();

    virtual void setLSData(int ls_idx, int vol_idx);
    virtual void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);
    virtual inline void setStencilWidth(int stencil_size)
    {
        d_update_weights = true;
        d_stencil_size = stencil_size;
    }

    inline int getStencilWidth()
    {
        return d_stencil_size;
    }

    virtual inline void setUseCentroids(bool use_centroids)
    {
        d_use_centroids = use_centroids;
        d_update_weights = true;
    }

    virtual void cacheData() = 0;

    virtual double reconstructOnIndex(IBTK::VectorNd x_loc,
                                      const SAMRAI::hier::Index<NDIM>& idx,
                                      const SAMRAI::pdat::CellData<NDIM, double>& Q_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch) = 0;

protected:
    int d_stencil_size = 2;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    std::vector<std::map<PatchIndexPair, Eigen::FullPivHouseholderQR<IBTK::MatrixXd>>> d_qr_matrix_vec;
    bool d_update_weights = true;
    bool d_use_centroids = true;
    int d_vol_idx = IBTK::invalid_index, d_ls_idx = IBTK::invalid_index;
};
} // namespace LS
#endif
