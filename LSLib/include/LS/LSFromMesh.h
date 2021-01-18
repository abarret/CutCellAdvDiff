#ifndef included_LS_LSFromMesh
#define included_LS_LSFromMesh

#include "ibtk/config.h"

#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IndexUtilities.h"

#include "LS/LSFindCellVolume.h"

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

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

namespace LS
{
class LSFromMesh : public LSFindCellVolume
{
public:
    LSFromMesh(std::string object_name,
               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
               libMesh::MeshBase* mesh,
               IBTK::FEDataManager* fe_data_manager,
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

private:
    bool findIntersection(libMesh::Point& p, libMesh::Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q);

    libMesh::MeshBase* d_mesh;
    IBTK::FEDataManager* d_fe_data_manager;
    bool d_use_inside = true;

    static const double s_eps;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_sgn_var;
    int d_sgn_idx = IBTK::invalid_index;
};
} // namespace LS

#endif
