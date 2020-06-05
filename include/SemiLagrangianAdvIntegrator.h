#ifndef included_SemiLagrangianAdvIntegrator
#define included_SemiLagrangianAdvIntegrator

#include "IBAMR_config.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/LSInitStrategy.h"

#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/PoissonSolver.h"

#include "LSCutCellLaplaceOperator.h"
#include "LSFindCellVolume.h"

#include "utility_functions.h"

namespace IBAMR
{
class SemiLagrangianAdvIntegrator : public AdvDiffHierarchyIntegrator
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

    void registerLevelSetVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var);

    void registerLevelSetVelocity(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
                                  SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double>> u_var);

    void registerLevelSetResetFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
                                       SAMRAI::tbox::Pointer<LSInitStrategy> ls_strategy);

    void registerLevelSetFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
                                  SAMRAI::tbox::Pointer<IBTK::CartGridFunction> ls_fcn);

    void restrictToLevelSet(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var);

    void useLevelSetFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
                             bool use_ls_function);

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

    void initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM>> old_level,
                                        bool allocate_data) override;

    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                          int ln,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too) override;

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
    void regridHierarchyEndSpecialized() override;
    void resetTimeDependentHierarchyDataSpecialized(double new_time) override;

private:
    void advectionUpdate(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                         double current_time,
                         double new_time);

    void diffusionUpdate(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                         int ls_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                         int vol_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vol_var,
                         int area_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> area_var,
                         double current_time,
                         double new_time);

    void integratePaths(int path_idx, int u_idx, double dt);
    void integratePaths(int path_idx, int u_idx, int vol_idx, int ls_idx, double dt);

    void evaluateMappingOnHierarchy(int xstar_idx,
                                    int Q_cur_idx,
                                    int vol_cur_idx,
                                    int Q_new_idx,
                                    int vol_new_idx,
                                    int ls_idx,
                                    int order);

    void evaluateMappingOnHierarchy(int xstar_idx, int Q_cur_idx, int Q_new_idx, int order);

    double sumOverZSplines(const IBTK::VectorNd& x_loc,
                           const SAMRAI::pdat::CellIndex<NDIM>& idx,
                           const SAMRAI::pdat::CellData<NDIM, double>& Q_data,
                           const int order);

    bool indexWithinWidth(int stencil_width,
                          const CellIndex<NDIM>& idx,
                          const SAMRAI::pdat::CellData<NDIM, double>& vol_data);

    double leastSquaresReconstruction(IBTK::VectorNd x_loc,
                                      const CellIndex<NDIM>& idx,
                                      const CellData<NDIM, double>& Q_data,
                                      const CellData<NDIM, double>& vol_data,
                                      const NodeData<NDIM, double>& ls_data,
                                      const Pointer<Patch<NDIM>>& patch);

    VectorNd
    findVelocity(const CellIndex<NDIM>& idx, const SAMRAI::pdat::SideData<NDIM, double>& u_data, const VectorNd& x_pt);

    double evaluateZSpline(const IBTK::VectorNd x, const int order);

    int getSplineWidth(int order);
    double ZSpline(double x, int order);

    double weight(double r);

    // patch data for particle trajectories
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_path_var;
    int d_path_idx = IBTK::invalid_index;
    int d_Q_scratch_idx = IBTK::invalid_index;

    SAMRAI::hier::ComponentSelector d_adv_data;

    // Level set information
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>> d_ls_node_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>> d_ls_cell_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>> d_vol_vars, d_area_vars;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>,
             SAMRAI::tbox::Pointer<IBTK::CartGridFunction>>
        d_ls_fcn_map;
    SAMRAI::tbox::Pointer<LSFindCellVolume> d_vol_fcn;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>, bool> d_ls_use_fcn;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>>
        d_Q_ls_map;

    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>, SAMRAI::tbox::Pointer<LSInitStrategy>>
        d_ls_strategy_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double>>>
        d_ls_u_map;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_u_s_var;

    SAMRAI::hier::ComponentSelector d_ls_data;

    static int GHOST_CELL_WIDTH;

    double d_min_ls_refine_factor = std::numeric_limits<double>::quiet_NaN();
    double d_max_ls_refine_factor = std::numeric_limits<double>::quiet_NaN();
    LeastSquaresOrder d_least_squares_reconstruction_order = UNKNOWN_ORDER;
    bool d_use_strang_splitting = false;

}; // Class SemiLagrangianAdvIntegrator
} // Namespace IBAMR

#include "private/SemiLagrangianAdvIntegrator_inc.h"

#endif
