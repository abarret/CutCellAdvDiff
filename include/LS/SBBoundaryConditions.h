#ifndef included_SBBoundaryConditions
#define included_SBBoundaryConditions

#include "ibtk/FEDataManager.h"

#include "LS/LSCutCellBoundaryConditions.h"
#include "LS/RBFReconstructCache.h"
#include "LS/SBSurfaceFluidCouplingManager.h"

#include "libmesh/equation_systems.h"
#include <libmesh/mesh.h>

namespace LS
{
class SBBoundaryConditions : public LSCutCellBoundaryConditions
{
public:
    using ReactionFcn =
        std::function<double(double, const std::vector<double>&, const std::vector<double>&, double, void*)>;

    SBBoundaryConditions(const std::string& object_name,
                         const std::string& fl_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         libMesh::Mesh* mesh,
                         const std::shared_ptr<SBSurfaceFluidCouplingManager>& sb_data_manager,
                         const std::shared_ptr<CutCellMeshMapping>& cut_cell_mesh_mapping);

    ~SBBoundaryConditions() = default;

    /*!
     * \brief Deleted default constructor.
     */
    SBBoundaryConditions() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    SBBoundaryConditions(const SBBoundaryConditions& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    SBBoundaryConditions& operator=(const SBBoundaryConditions& that) = delete;

    void setFluidContext(SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx);

    void allocateOperatorState(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                               double time) override;
    void deallocateOperatorState(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double time) override;

    void applyBoundaryCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                int Q_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> R_var,
                                int R_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                double time) override;

private:

    double d_D_coef = std::numeric_limits<double>::quiet_NaN();

    libMesh::Mesh* d_mesh = nullptr;
    std::shared_ptr<SBSurfaceFluidCouplingManager> d_sb_data_manager;
    std::shared_ptr<CutCellMeshMapping> d_cut_cell_mapping;

    std::string d_fl_name;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;
};

} // namespace LS
#endif
