#ifndef included_SBBoundaryConditions
#define included_SBBoundaryConditions

#include "ibtk/FECache.h"
#include "ibtk/FEDataManager.h"

#include "LS/LSCutCellBoundaryConditions.h"
#include "LS/RBFReconstructCache.h"

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
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         libMesh::Mesh* mesh,
                         IBTK::FEDataManager* fe_manager,
                         bool use_old_soln = true);

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

    void registerFluidFluidInteraction(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> fl_var);
    void setFluidContext(SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx);

    void registerFluidSurfaceInteraction(const std::string& surface_name);

    void allocateOperatorState(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                               double time) override;
    void deallocateOperatorState(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double time) override;

    inline void setReactionFunction(ReactionFcn a_fcn, ReactionFcn g_fcn, void* ctx)
    {
        d_a_fcn = a_fcn;
        d_g_fcn = g_fcn;
        d_fcn_ctx = ctx;
    }

    void applyBoundaryCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                int Q_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> R_var,
                                int R_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                double time) override;

private:
    void interpolateToBoundary(int Q_idx,
                               const std::string& Q_sys_name,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                               double current_time);

    /*!
     * Find an intersection between the element elem and the side defined by the point r and the search direction and
     * magnitude q.
     */
    bool findIntersection(libMesh::Point& p, libMesh::Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q);

    void clearCache();

    void cacheBoundaryData();

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    std::vector<std::map<LS::PatchIndexPair, std::vector<std::unique_ptr<libMesh::Node>>>> d_node_idx_map_vec;
    std::vector<std::map<LS::PatchIndexPair, std::vector<std::pair<libMesh::Elem*, std::unique_ptr<libMesh::Elem>>>>>
        d_elem_idx_map_vec;
    std::vector<std::map<LS::PatchIndexPair, std::vector<std::vector<libMesh::Point>>>> d_intersection_idx_map_vec;

    double d_D_coef = std::numeric_limits<double>::quiet_NaN();

    libMesh::Mesh* d_mesh = nullptr;
    IBTK::FEDataManager* d_fe_data_manager = nullptr;
    std::vector<std::string> d_sf_names;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>> d_fl_vars;
    std::vector<std::string> d_fl_names;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;
    std::string d_sys_name;
    ReactionFcn d_a_fcn, d_g_fcn;
    void* d_fcn_ctx = nullptr;

    bool d_perturb_nodes = true;
    bool d_use_old_soln = true;

    RBFReconstructCache d_rbf_reconstruct;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_scr_var;
    int d_scr_idx = IBTK::invalid_index;
};

} // namespace LS
#endif
