#ifndef included_SBIntegrator
#define included_SBIntegrator

#include "ibtk/FEDataManager.h"

#include "LS/CutCellMeshMapping.h"
#include "LS/RBFReconstructCache.h"
#include "LS/SBSurfaceFluidCouplingManager.h"

#include "CellData.h"
#include "NodeData.h"

#include "libmesh/equation_systems.h"
#include <libmesh/mesh.h>

namespace LS
{
class SBIntegrator
{
public:
    using ReactionFcn =
        std::function<double(double, const std::vector<double>&, const std::vector<double>&, double, void*)>;

    SBIntegrator(std::string object_name,
                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                 const std::shared_ptr<SBSurfaceFluidCouplingManager>& sb_sf_fl_coupling_manager,
                 libMesh::Mesh* mesh);

    /*!
     * \brief Deleted default constructor.
     */
    SBIntegrator() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    SBIntegrator(const SBIntegrator& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    SBIntegrator& operator=(const SBIntegrator& that) = delete;

    void setLSData(int ls_idx, int vol_idx, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    void
    integrateHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx, double current_time, double new_time);

    void beginTimestepping(double current_time, double new_time);
    void endTimestepping(double current_time, double new_time);

    void interpolateToBoundary(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> fl_var,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
                               double time);

private:
    std::string d_object_name;

    std::shared_ptr<SBSurfaceFluidCouplingManager> d_sb_data_manager = nullptr;

    libMesh::Mesh* d_mesh = nullptr;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    int d_vol_idx = IBTK::invalid_index, d_ls_idx = IBTK::invalid_index;
};

} // namespace LS
#endif
