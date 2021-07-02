#ifndef included_CCAD_LSFromMesh
#define included_CCAD_LSFromMesh

#include "ibtk/config.h"

#include "CCAD/CutCellMeshMapping.h"
#include "CCAD/FEMeshPartitioner.h"
#include "CCAD/LSFindCellVolume.h"

#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/app_namespaces.h"

#include "Box.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Pointer.h"

#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/vector_value.h"

#include <functional>

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

namespace CCAD
{
class LSFromMesh : public LSFindCellVolume
{
public:
    LSFromMesh(std::string object_name,
               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
               const SAMRAI::tbox::Pointer<CutCellMeshMapping>& cut_cell_mesh_mapping,
               bool use_inside = true);

    virtual ~LSFromMesh() = default;

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

    inline void registerNormalReverseDomainId(unsigned int bdry_id, unsigned int part = 0)
    {
        d_norm_reverse_domain_ids[part].insert(bdry_id);
    }

    inline void registerNormalReverseDomainId(const std::vector<unsigned int>& bdry_ids, unsigned int part = 0)
    {
        for (const auto& bdry_id : bdry_ids) registerNormalReverseDomainId(bdry_id, part);
    }

    inline void registerNormalReverseElemId(unsigned int bdry_id, unsigned int part = 0)
    {
        d_norm_reverse_elem_ids[part].insert(bdry_id);
    }

    inline void registerNormalReverseElemId(const std::vector<unsigned int>& bdry_ids, unsigned int part = 0)
    {
        for (const auto& bdry_id : bdry_ids) registerNormalReverseElemId(bdry_id, part);
    }

    using BdryFcn = std::function<void(const IBTK::VectorNd&, double&)>;

    inline void registerBdryFcn(BdryFcn fcn)
    {
        d_bdry_fcn = fcn;
    }

private:
    void commonConstructor();
    using Simplex = std::array<std::pair<IBTK::VectorNd, double>, NDIM + 1>;
    using LDSimplex = std::array<std::pair<IBTK::VectorNd, double>, NDIM>;
    using PolytopePt = std::tuple<IBTK::VectorNd, int, int>;
    void findVolume(const double* const xlow,
                    const double* const dx,
                    const SAMRAI::hier::Index<NDIM>& patch_lower,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double>> phi_data,
                    const SAMRAI::pdat::CellIndex<NDIM>& idx,
                    double& volume);
    double findVolume(const std::vector<Simplex>& simplices);

    void updateLSAwayFromInterface(int phi_idx);
    // This does a flood filling algorithm for d_sgn_idx.
    // We assume that any value less than eps on the given level is correctly set.
    // NOTE: eps must be positive.
    void floodFillForLS(int ln, double eps);

    bool d_use_inside = true;

    SAMRAI::tbox::Pointer<CutCellMeshMapping> d_cut_cell_mesh_mapping;

    std::vector<std::set<unsigned int>> d_norm_reverse_domain_ids, d_norm_reverse_elem_ids;

    BdryFcn d_bdry_fcn;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_sgn_var;
    int d_sgn_idx = IBTK::invalid_index;
};
} // namespace CCAD

#endif
