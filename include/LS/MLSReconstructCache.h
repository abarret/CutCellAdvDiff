#ifndef included_MLSReconstructCache
#define included_MLSReconstructCache

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
 * \brief Class MLSReconstructCache caches the data necessary to form RBF reconstructions of data.
 */
class MLSReconstructCache : public ReconstructCache
{
public:
    MLSReconstructCache() = default;

    MLSReconstructCache(int ls_idx,
                        int vol_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                        bool use_centroids = true);

    ~MLSReconstructCache() = default;

    /*!
     * \brief Deleted copy constructor.
     */
    MLSReconstructCache(const MLSReconstructCache& from) = delete;

    void cacheData() override;

    double reconstructOnIndex(IBTK::VectorNd x_loc,
                              const SAMRAI::hier::Index<NDIM>& idx,
                              const SAMRAI::pdat::CellData<NDIM, double>& Q_data,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch);
};
} // namespace LS
#endif
