#ifndef included_LS_SemiLagrangianAdvIntegrator
#define included_LS_SemiLagrangianAdvIntegrator

#include "ibamr/config.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/LSInitStrategy.h"

#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/PoissonSolver.h"

#include "LS/LSCutCellLaplaceOperator.h"
#include "LS/LSFindCellVolume.h"
#include "LS/MLSReconstructCache.h"
#include "LS/RBFReconstructCache.h"
#include "LS/SBIntegrator.h"
#include "LS/utility_functions.h"

namespace LS
{
class SemiLagrangianAdvIntegrator : public IBAMR::AdvDiffHierarchyIntegrator
{
public:
    SemiLagrangianAdvIntegrator(const std::string& object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                bool register_for_restart = false);

    ~SemiLagrangianAdvIntegrator() = default;

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator. Can optionally turn off outputting the quantity.
     *
     * Data management for the registered quantity will be handled by the
     * hierarchy integrator.
     */
    void registerTransportedQuantity(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                     const bool Q_output = true) override;

    void registerLevelSetVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    void registerLevelSetVolFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                     SAMRAI::tbox::Pointer<LSFindCellVolume> vol_fcn);

    void registerLevelSetResetFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                       SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_strategy);

    void setFEDataManagerNeedsInitialization(IBTK::FEDataManager* fe_data_manager);

    void registerLevelSetSBDataManager(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                       std::shared_ptr<SBSurfaceFluidCouplingManager> sb_data_manager);

    void registerLevelSetCutCellMapping(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                        std::shared_ptr<CutCellMeshMapping> cut_cell_mesh_mapping);

    void restrictToLevelSet(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    void useLevelSetForTagging(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                               bool use_ls_for_tagging);

    void evolveLevelSet(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double>> u_var);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>
    getAreaVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>
    getVolumeVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>
    getLSCellVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    void registerSBIntegrator(SAMRAI::tbox::Pointer<SBIntegrator> sb_integrator,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg) override;

    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                          int ln,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too) override;

    void initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM>> old_level,
                                        bool allocate_data) override;

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int getNumberOfCycles() const override;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

protected:
    void initializeCompositeHierarchyDataSpecialized(double current_time, bool initial_time) override;
    void regridHierarchyBeginSpecialized() override;
    void resetTimeDependentHierarchyDataSpecialized(double new_time) override;
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> base_hierarchy,
                                           int coarsest_ln,
                                           int finest_ln);

    void addWorkloadEstimate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                             const int workload_data_idx) override;

    virtual void advectionUpdate(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                 double current_time,
                                 double new_time);

    virtual void diffusionUpdate(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                 int ls_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                 int vol_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vol_var,
                                 int area_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> area_var,
                                 int side_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> side_var,
                                 double current_time,
                                 double new_time);

    virtual void integratePaths(int path_idx, int u_new_idx, int u_half_idx, double dt);
    virtual void integratePaths(int path_idx, int u_new_idx, int u_half_idx, int vol_idx, int ls_idx, double dt);

    virtual void evaluateMappingOnHierarchy(int xstar_idx,
                                            int Q_cur_idx,
                                            int vol_cur_idx,
                                            int Q_new_idx,
                                            int vol_new_idx,
                                            int ls_idx,
                                            int order);

    SAMRAI::hier::ComponentSelector d_adv_data;

    // Level set information
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>> d_ls_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>> d_vol_vars, d_area_vars,
        d_vol_wgt_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>>> d_side_vars;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>, bool> d_ls_use_ls_for_tagging;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>>
        d_Q_ls_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>, SAMRAI::tbox::Pointer<LSFindCellVolume>>
        d_ls_vol_fcn_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double>>>
        d_ls_u_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>>
        d_ls_ls_cell_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>,
             SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy>>
        d_ls_strategy_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>,
             std::shared_ptr<SBSurfaceFluidCouplingManager>>
        d_ls_sb_data_manager_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>, std::shared_ptr<CutCellMeshMapping>>
        d_ls_cut_cell_mapping_map;
    std::vector<IBTK::FEDataManager*> d_fe_data_managers;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_u_s_var;

    // patch data for particle trajectories
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_path_var;
    int d_path_idx = IBTK::invalid_index;
    // Scratch data for when we need more ghost cells.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_Q_big_scr_var;
    int d_Q_big_scr_idx = IBTK::invalid_index;

    SAMRAI::hier::ComponentSelector d_ls_data;

    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>> d_sol_ls_vecs, d_rhs_ls_vecs;

    static int GHOST_CELL_WIDTH;

    double d_min_ls_refine_factor = std::numeric_limits<double>::quiet_NaN();
    double d_max_ls_refine_factor = std::numeric_limits<double>::quiet_NaN();
    LeastSquaresOrder d_least_squares_reconstruction_order = LeastSquaresOrder::UNKNOWN_ORDER;
    AdvectionTimeIntegrationMethod d_adv_ts_type = AdvectionTimeIntegrationMethod::UNKNOWN_METHOD;
    DiffusionTimeIntegrationMethod d_dif_ts_type = DiffusionTimeIntegrationMethod::UNKNOWN_METHOD;
    bool d_use_strang_splitting = false;

    bool d_use_rbfs = false;
    int d_rbf_stencil_size = 2;

    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double>> d_hier_fc_data_ops;

    std::map<SAMRAI::tbox::Pointer<SBIntegrator>, SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>>
        d_sb_integrator_ls_map;

    SAMRAI::tbox::Pointer<ReconstructCache> d_reconstruction_cache;

private:
    void evaluateMappingOnHierarchy(int xstar_idx, int Q_cur_idx, int Q_new_idx, int order);

    double sumOverZSplines(const IBTK::VectorNd& x_loc,
                           const SAMRAI::pdat::CellIndex<NDIM>& idx,
                           const SAMRAI::pdat::CellData<NDIM, double>& Q_data,
                           const int order);

    bool indexWithinWidth(int stencil_width,
                          const SAMRAI::pdat::CellIndex<NDIM>& idx,
                          const SAMRAI::pdat::CellData<NDIM, double>& vol_data);

    double radialBasisFunctionReconstruction(IBTK::VectorNd x_loc,
                                             const SAMRAI::pdat::CellIndex<NDIM>& idx,
                                             const SAMRAI::pdat::CellData<NDIM, double>& Q_data,
                                             const SAMRAI::pdat::CellData<NDIM, double>& vol_data,
                                             const SAMRAI::pdat::NodeData<NDIM, double>& ls_data,
                                             const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch);

    double leastSquaresReconstruction(IBTK::VectorNd x_loc,
                                      const SAMRAI::pdat::CellIndex<NDIM>& idx,
                                      const SAMRAI::pdat::CellData<NDIM, double>& Q_data,
                                      const SAMRAI::pdat::CellData<NDIM, double>& vol_data,
                                      const SAMRAI::pdat::NodeData<NDIM, double>& ls_data,
                                      const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch);

    double evaluateZSpline(const IBTK::VectorNd x, const int order);

    int getSplineWidth(int order);
    double ZSpline(double x, int order);

    double weight(double r);

    RBFPolyOrder d_rbf_poly_order = RBFPolyOrder::UNKNOWN_ORDER;
}; // Class SemiLagrangianAdvIntegrator
} // Namespace LS

#include "LS/private/SemiLagrangianAdvIntegrator_inc.h"

#endif
