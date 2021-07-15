#ifndef included_CCAD_LSAdvDiffIntegrator
#define included_CCAD_LSAdvDiffIntegrator

#include "ibamr/config.h"

#include "CCAD/AdvectiveReconstructionOperator.h"
#include "CCAD/LSCutCellLaplaceOperator.h"
#include "CCAD/LSFindCellVolume.h"
#include "CCAD/MLSReconstructCache.h"
#include "CCAD/RBFReconstructCache.h"
#include "CCAD/SBIntegrator.h"
#include "CCAD/VolumeBoundaryMeshMapping.h"
#include "CCAD/ls_utilities.h"
#include "CCAD/reconstructions.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/LSInitStrategy.h"

#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/PoissonSolver.h"

namespace CCAD
{
class LSAdvDiffIntegrator : public IBAMR::AdvDiffHierarchyIntegrator
{
public:
    LSAdvDiffIntegrator(const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        bool register_for_restart = true);

    ~LSAdvDiffIntegrator() = default;

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator. Can optionally turn off outputting the quantity.
     *
     * Data management for the registered quantity will be handled by the
     * hierarchy integrator.
     */
    void registerTransportedQuantity(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                     const bool Q_output = true) override;

    virtual void registerLevelSetVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    virtual void registerLevelSetVolFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                             SAMRAI::tbox::Pointer<LSFindCellVolume> vol_fcn);

    virtual void registerLevelSetResetFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                               SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_strategy);

    virtual void restrictToLevelSet(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    virtual void useLevelSetForTagging(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                       bool use_ls_for_tagging);

    virtual void evolveLevelSet(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double>> u_var);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>
    getAreaVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>
    getVolumeVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>
    getLSCellVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

    void registerReconstructionCache(SAMRAI::tbox::Pointer<CCAD::ReconstructCache> reconstruct_cache);

    void registerAdvectionReconstruction(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                         std::shared_ptr<AdvectiveReconstructionOperator> reconstruct_op);

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
    void resetTimeDependentHierarchyDataSpecialized(double new_time) override;
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> base_hierarchy,
                                           int coarsest_ln,
                                           int finest_ln) override;

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

    virtual void evaluateMappingOnHierarchy(int xstar_idx, int Q_cur_idx, int Q_new_idx, int order);

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

    // Advection reconstruction information
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>,
             std::shared_ptr<AdvectiveReconstructionOperator>>
        d_Q_adv_reconstruct_map;
    AdvReconstructType d_default_adv_reconstruct_type = AdvReconstructType::RBF;
    std::shared_ptr<AdvectiveReconstructionOperator> d_default_adv_reconstruct_op;

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
    Reconstruct::LeastSquaresOrder d_least_squares_reconstruction_order = Reconstruct::LeastSquaresOrder::UNKNOWN_ORDER;
    AdvectionTimeIntegrationMethod d_adv_ts_type = AdvectionTimeIntegrationMethod::UNKNOWN_METHOD;
    DiffusionTimeIntegrationMethod d_dif_ts_type = DiffusionTimeIntegrationMethod::UNKNOWN_METHOD;
    bool d_use_strang_splitting = false;

    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double>> d_hier_fc_data_ops;

    SAMRAI::tbox::Pointer<ReconstructCache> d_reconstruction_cache;

private:
    bool d_use_rbfs = false;
    unsigned int d_rbf_stencil_size = 2;
    Reconstruct::RBFPolyOrder d_rbf_poly_order = Reconstruct::RBFPolyOrder::UNKNOWN_ORDER;
}; // Class LSAdvDiffIntegrator
} // namespace CCAD

#endif
