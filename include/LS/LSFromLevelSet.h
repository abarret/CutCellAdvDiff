#ifndef included_LS_LSFromLevelSet
#define included_LS_LSFromLevelSet

#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_utilities.h"

#include "LS/LSFindCellVolume.h"

#include "Box.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "tbox/Pointer.h"

#include <array>
#include <vector>

namespace LS
{
class LSFromLevelSet : public LSFindCellVolume
{
public:
    LSFromLevelSet(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    ~LSFromLevelSet() = default;

    void registerLSFcn(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> ls_fcn);

    void updateVolumeAreaSideLS(int vol_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vol_var,
                                int area_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> area_var,
                                int side_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> side_var,
                                int phi_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> phi_var,
                                double data_time,
                                bool extended_box = false) override;

private:
    // Volume calculations
    using Simplex = std::array<std::pair<IBTK::VectorNd, double>, NDIM + 1>;
    using LDSimplex = std::array<std::pair<IBTK::VectorNd, double>, NDIM>;
    using PolytopePt = std::tuple<IBTK::VectorNd, int, int>;
    void findVolumeAndArea(const double* const xlow,
                           const double* const dx,
                           const SAMRAI::hier::Index<NDIM>& patch_lower,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double>> phi,
                           const SAMRAI::pdat::CellIndex<NDIM>& idx,
                           double& volume,
                           double& area);
    double findVolume(const std::vector<Simplex>& simplices);
    double findArea(const std::vector<Simplex>& simplices);

    static const double s_eps;

    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_ls_fcn;
};
} // namespace LS

#endif
