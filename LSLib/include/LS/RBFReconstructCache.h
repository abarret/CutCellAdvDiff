#ifndef included_RBFReconstructCache
#define included_RBFReconstructCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"

#include "LS/ReconstructCache.h"
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
 * \brief Class RBFReconstructCache caches the data necessary to form RBF reconstructions of data.
 */
class RBFReconstructCache : public ReconstructCache
{
public:
    RBFReconstructCache() = default;

    RBFReconstructCache(int ls_idx,
                        int vol_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                        bool use_centroids = true);

    ~RBFReconstructCache() = default;

    /*!
     * \brief Deleted copy constructor.
     */
    RBFReconstructCache(const RBFReconstructCache& from) = delete;

    void cacheData() override;

    double reconstructOnIndex(IBTK::VectorNd x_loc,
                              const SAMRAI::hier::Index<NDIM>& idx,
                              const SAMRAI::pdat::CellData<NDIM, double>& Q_data,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch) override;
};
} // namespace LS
#endif
