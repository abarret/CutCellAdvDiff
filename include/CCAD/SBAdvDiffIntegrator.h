#ifndef included_CCAD_SBAdvDiffIntegrator
#define included_CCAD_SBAdvDiffIntegrator

#include "ibamr/config.h"

#include "CCAD/AdvectiveReconstructionOperator.h"
#include "CCAD/LSAdvDiffIntegrator.h"
#include "CCAD/LSCutCellLaplaceOperator.h"
#include "CCAD/LSFindCellVolume.h"
#include "CCAD/MLSReconstructCache.h"
#include "CCAD/RBFReconstructCache.h"
#include "CCAD/SBIntegrator.h"
#include "CCAD/VolumeBoundaryMeshMapping.h"
#include "CCAD/ls_utilities.h"
#include "CCAD/reconstructions.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/LSInitStrategy.h"

#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/PoissonSolver.h"

namespace CCAD
{
void callback_fcn(double current_time, double new_time, int cycle_num, void* ctx);

class SBAdvDiffIntegrator : public LSAdvDiffIntegrator
{
public:
    SBAdvDiffIntegrator(const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        SAMRAI::tbox::Pointer<IBAMR::IBHierarchyIntegrator> ib_integrator,
                        bool register_for_restart = true);

    ~SBAdvDiffIntegrator() = default;

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator. Can optionally turn off outputting the quantity.
     *
     * Data management for the registered quantity will be handled by the
     * hierarchy integrator.
     */

    static void callbackIntegrateHierarchy(double current_time, double neW_time, int cycle_num, void* ctx);

    void registerVolumeBoundaryMeshMapping(const std::shared_ptr<VolumeBoundaryMeshMapping>& vol_bdry_mesh_mapping);

    void registerLevelSetSBDataManager(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                       std::shared_ptr<SBSurfaceFluidCouplingManager> sb_data_manager);

    void registerLevelSetCutCellMapping(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var,
                                        std::shared_ptr<CutCellMeshMapping> cut_cell_mesh_mapping);

    void registerSBIntegrator(SAMRAI::tbox::Pointer<SBIntegrator> sb_integrator,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> ls_var);

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

protected:
    void initializeCompositeHierarchyDataSpecialized(double current_time, bool initial_time) override;
    void regridHierarchyEndSpecialized() override;

    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>,
             std::shared_ptr<SBSurfaceFluidCouplingManager>>
        d_ls_sb_data_manager_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>, std::shared_ptr<CutCellMeshMapping>>
        d_ls_cut_cell_mapping_map;
    std::shared_ptr<VolumeBoundaryMeshMapping> d_vol_bdry_mesh_mapping;
    std::map<SAMRAI::tbox::Pointer<SBIntegrator>, SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>>
        d_sb_integrator_ls_map;

private:
    bool d_use_rbfs = false;
    unsigned int d_rbf_stencil_size = 2;
    Reconstruct::RBFPolyOrder d_rbf_poly_order = Reconstruct::RBFPolyOrder::UNKNOWN_ORDER;
}; // Class SBAdvDiffIntegrator

} // namespace CCAD

#endif
