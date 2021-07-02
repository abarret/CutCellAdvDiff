#include "ibamr/config.h"

#include "CCAD/LSAdvDiffIntegrator.h"
#include "CCAD/LSCutCellLaplaceOperator.h"
#include "CCAD/LSFindCellVolume.h"
#include "CCAD/LSFromLevelSet.h"
#include "CCAD/app_namespaces.h"
#include "CCAD/ls_functions.h"

#include <ibamr/RelaxationLSMethod.h>

#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/muParserCartGridFunction.h"
#include "ibtk/muParserRobinBcCoefs.h"
#include <ibtk/AppInitializer.h>

#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "LoadBalancer.h"
#include "SAMRAI_config.h"
#include "StandardTagAndInitialize.h"
#include "tbox/Pointer.h"

#include <petscsys.h>

#include <memory>
#include <utility>
// Local includes
#include "InsideBoundaryConditions.h"
#include "InsideLSFcn.h"
#include "OutsideBoundaryConditions.h"
#include "OutsideLSFcn.h"
#include "QFcn.h"

static double a = std::numeric_limits<double>::signaling_NaN();
static double b = std::numeric_limits<double>::signaling_NaN();

void postprocess_data(Pointer<PatchHierarchy<NDIM>> hierarchy,
                      Pointer<LSAdvDiffIntegrator> integrator,
                      Pointer<CellVariable<NDIM, double>> Q_in_var,
                      Pointer<CellVariable<NDIM, double>> Q_out_var,
                      int iteration_num,
                      double loop_time,
                      const std::string& dirname);

void output_to_file(const int Q_idx,
                    const int area_idx,
                    const int vol_idx,
                    const int ls_interior_idx,
                    const std::string& file_name,
                    const double loop_time,
                    Pointer<PatchHierarchy<NDIM>> hierarchy);
void
outputBdryInfo(const int Q_idx,
               const int Q_scr_idx,
               const int ls_interior_idx,
               const int vol_idx,
               const int area_idx,
               const double current_time,
               const int iteration_num,
               Pointer<PatchHierarchy<NDIM>> hierarchy,
               std::string base_name)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(Q_scr_idx)) level->allocatePatchData(Q_scr_idx);
    }

    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comps(1);
    ghost_cell_comps[0] = ITC(Q_scr_idx, Q_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr);
    HierarchyGhostCellInterpolation hier_ghost_cells;
    hier_ghost_cells.initializeOperatorState(ghost_cell_comps, hierarchy, 0, hierarchy->getFinestLevelNumber());
    hier_ghost_cells.fillData(current_time);

    output_to_file(Q_scr_idx,
                   area_idx,
                   vol_idx,
                   ls_interior_idx,
                   base_name + std::to_string(iteration_num),
                   current_time,
                   hierarchy);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(Q_scr_idx);
    }
}

void
calculateTotalAmounts(Pointer<PatchHierarchy<NDIM>> hierarchy,
                      const double time,
                      const int Q_in_idx,
                      const int Q_out_idx,
                      const int vol_in_idx,
                      const int vol_out_idx)
{
    std::ofstream amounts_stream;
    if (SAMRAI_MPI::getRank() == 0)
    {
        if (time == 0.0)
            amounts_stream.open("total_amount", std::ofstream::out);
        else
            amounts_stream.open("total_amount", std::ofstream::app);
    }
    Pointer<HierarchyMathOps> hier_math_ops = new HierarchyMathOps("HierarchyMathOps", hierarchy);
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
    double tot_out = 0.0, tot_in = 0.0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();

            Pointer<CellData<NDIM, double>> Q_in_data = patch->getPatchData(Q_in_idx);
            Pointer<CellData<NDIM, double>> Q_out_data = patch->getPatchData(Q_out_idx);
            Pointer<CellData<NDIM, double>> vol_in_data = patch->getPatchData(vol_in_idx);
            Pointer<CellData<NDIM, double>> vol_out_data = patch->getPatchData(vol_out_idx);
            Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);

            for (CellIterator<NDIM> ci(box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                tot_out += (*Q_out_data)(idx) * (*wgt_data)(idx) * (*vol_out_data)(idx);
                tot_in += (*Q_in_data)(idx) * (*wgt_data)(idx) * (*vol_in_data)(idx);
            }
        }
    }
    tot_out = SAMRAI_MPI::sumReduction(tot_out);
    tot_in = SAMRAI_MPI::sumReduction(tot_in);
    pout << std::setprecision(10) << "Total interior at time t = " << time << " is " << tot_in << "\n";
    pout << std::setprecision(10) << "Total exterior at time t = " << time << " is " << tot_out << "\n";
    pout << std::setprecision(10) << "Total amount at time t =   " << time << " is " << tot_out + tot_in << "\n";
    if (SAMRAI_MPI::getRank() == 0)
    {
        amounts_stream << time << " " << std::setprecision(10) << tot_in << " " << tot_out << " " << tot_out + tot_in
                       << "\n";
        amounts_stream.close();
    }
}

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int dump_postproc_interval = app_initializer->getPostProcessingDataDumpInterval();
        const std::string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<LSAdvDiffIntegrator> time_integrator = new LSAdvDiffIntegrator(
            "LSAdvDiffIntegrator", app_initializer->getComponentDatabase("LSAdvDiffIntegrator"), false);

        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Setup the advection velocity.
        Pointer<FaceVariable<NDIM, double>> u_var = new FaceVariable<NDIM, double>("u");
        Pointer<CartGridFunction> u_fcn;
        Pointer<CartGridFunction> u_adv_fcn;
        time_integrator->registerAdvectionVelocity(u_var);
        time_integrator->setAdvectionVelocityIsDivergenceFree(u_var, true);
        u_fcn = new muParserCartGridFunction(
            "UFunction", app_initializer->getComponentDatabase("UFunction"), grid_geometry);
        time_integrator->setAdvectionVelocityFunction(u_var, u_fcn);

        // Setup the level set function
        Pointer<NodeVariable<NDIM, double>> ls_in_var = new NodeVariable<NDIM, double>("LS_In");
        time_integrator->registerLevelSetVariable(ls_in_var);
        Pointer<InsideLSFcn> ls_in_fcn =
            new InsideLSFcn("InsideLSFcn", app_initializer->getComponentDatabase("InsideLSFcn"));
        Pointer<LSFromLevelSet> vol_in_fcn = new LSFromLevelSet("VolInFcn", patch_hierarchy);
        vol_in_fcn->registerLSFcn(ls_in_fcn);
        time_integrator->registerLevelSetVolFunction(ls_in_var, vol_in_fcn);

        // Setup second level set
        Pointer<NodeVariable<NDIM, double>> ls_out_var = new NodeVariable<NDIM, double>("LS_Out");
        time_integrator->registerLevelSetVariable(ls_out_var);
        Pointer<OutsideLSFcn> ls_out_fcn =
            new OutsideLSFcn("LS_OUT_FCN", time_integrator, ls_in_var, app_initializer->getComponentDatabase("LS_Out"));
        Pointer<LSFromLevelSet> vol_out_fcn = new LSFromLevelSet("VolOutFcn", patch_hierarchy);
        vol_out_fcn->registerLSFcn(ls_out_fcn);
        time_integrator->registerLevelSetVolFunction(ls_out_var, vol_out_fcn);
        time_integrator->useLevelSetForTagging(ls_out_var, input_db->getBool("USE_OUT_LS_FOR_TAGGING"));

        a = input_db->getDouble("A");
        b = input_db->getDouble("B");

        // Setup advected quantity
        Pointer<CellVariable<NDIM, double>> Q_in_var = new CellVariable<NDIM, double>("Q_in");
        Pointer<CartGridFunction> Q_in_init = new QFcn("QInit", app_initializer->getComponentDatabase("QFcnIn"));
        const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;
        std::vector<RobinBcCoefStrategy<NDIM>*> Q_in_bcs(1);
        if (!periodic_domain)
        {
            const std::string Q_bcs_name = "Q_bcs";
            Q_in_bcs[0] =
                new muParserRobinBcCoefs(Q_bcs_name, app_initializer->getComponentDatabase(Q_bcs_name), grid_geometry);
        }

        time_integrator->registerTransportedQuantity(Q_in_var);
        time_integrator->setAdvectionVelocity(Q_in_var, u_var);
        time_integrator->setInitialConditions(Q_in_var, Q_in_init);
        time_integrator->setPhysicalBcCoef(Q_in_var, Q_in_bcs[0]);
        time_integrator->setDiffusionCoefficient(Q_in_var, input_db->getDoubleWithDefault("D_coef", 0.0));
        time_integrator->restrictToLevelSet(Q_in_var, ls_in_var);

        // Setup advected quantity
        Pointer<CellVariable<NDIM, double>> Q_out_var = new CellVariable<NDIM, double>("Q_out");
        Pointer<QInitial> Q_out_init =
            new QInitial("QInit", grid_geometry, app_initializer->getComponentDatabase("QInitial_out"));
        std::vector<RobinBcCoefStrategy<NDIM>*> Q_out_bcs(1);
        if (!periodic_domain)
        {
            const std::string Q_bcs_name = "Q_bcs";
            Q_out_bcs[0] =
                new muParserRobinBcCoefs(Q_bcs_name, app_initializer->getComponentDatabase(Q_bcs_name), grid_geometry);
        }

        time_integrator->registerTransportedQuantity(Q_out_var);
        time_integrator->setAdvectionVelocity(Q_out_var, u_var);
        time_integrator->setInitialConditions(Q_out_var, Q_out_init);
        time_integrator->setPhysicalBcCoef(Q_out_var, Q_out_bcs[0]);
        time_integrator->setDiffusionCoefficient(Q_out_var, input_db->getDoubleWithDefault("D_coef", 0.0));
        time_integrator->restrictToLevelSet(Q_out_var, ls_out_var);

        // Set up diffusion operators
        Pointer<LSCutCellLaplaceOperator> rhs_in_oper = new LSCutCellLaplaceOperator(
            "LSCutCellRHSInOperator", app_initializer->getComponentDatabase("LSCutCellOperator"), false);
        Pointer<LSCutCellLaplaceOperator> sol_in_oper = new LSCutCellLaplaceOperator(
            "LSCutCellInOperator", app_initializer->getComponentDatabase("LSCutCellOperator"), false);
        time_integrator->setHelmholtzRHSOperator(Q_in_var, rhs_in_oper);
        Pointer<PETScKrylovPoissonSolver> Q_in_helmholtz_solver = new PETScKrylovPoissonSolver(
            "PoissonSolver", app_initializer->getComponentDatabase("PoissonSolver"), "poisson_solve_");
        Q_in_helmholtz_solver->setOperator(sol_in_oper);
        time_integrator->setHelmholtzSolver(Q_in_var, Q_in_helmholtz_solver);
        Pointer<InsideBoundaryConditions> in_bdry_oper = new InsideBoundaryConditions(
            "InsideBdryOper", app_initializer->getComponentDatabase("BdryConds"), Q_out_var, time_integrator);
        rhs_in_oper->setBoundaryConditionOperator(in_bdry_oper);
        sol_in_oper->setBoundaryConditionOperator(in_bdry_oper);

        Pointer<LSCutCellLaplaceOperator> rhs_out_oper = new LSCutCellLaplaceOperator(
            "LSCutCellRHSOutOperator", app_initializer->getComponentDatabase("LSCutCellOperator"), false);
        Pointer<LSCutCellLaplaceOperator> sol_out_oper = new LSCutCellLaplaceOperator(
            "LSCutCellOutOperator", app_initializer->getComponentDatabase("LSCutCellOperator"), false);
        time_integrator->setHelmholtzRHSOperator(Q_out_var, rhs_out_oper);
        Pointer<PETScKrylovPoissonSolver> Q_out_helmholtz_solver = new PETScKrylovPoissonSolver(
            "PoissonSolver", app_initializer->getComponentDatabase("PoissonSolver"), "poisson_solve_");
        Q_out_helmholtz_solver->setOperator(sol_out_oper);
        time_integrator->setHelmholtzSolver(Q_out_var, Q_out_helmholtz_solver);
        Pointer<OutsideBoundaryConditions> out_bdry_oper = new OutsideBoundaryConditions(
            "OutsideBdryOper", app_initializer->getComponentDatabase("BdryConds"), Q_in_var, time_integrator);
        rhs_out_oper->setBoundaryConditionOperator(out_bdry_oper);
        sol_out_oper->setBoundaryConditionOperator(out_bdry_oper);

        // Set up visualization plot file writer.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Register a drawing variable with the data writer
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<CellVariable<NDIM, double>> u_draw_var = new CellVariable<NDIM, double>("U", NDIM);
        const int u_draw_idx = var_db->registerVariableAndContext(u_draw_var, var_db->getContext("Draw"));
        visit_data_writer->registerPlotQuantity("U", "VECTOR", u_draw_idx);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        bool output_bdry_info = input_db->getBool("OUTPUT_BDRY_INFO");
        const int Q_out_idx = var_db->mapVariableAndContextToIndex(Q_out_var, time_integrator->getCurrentContext());
        const int Q_out_scr_idx =
            var_db->registerVariableAndContext(Q_out_var, var_db->getContext("SCRATCH"), IntVector<NDIM>(2));
        const int Q_in_idx = var_db->mapVariableAndContextToIndex(Q_in_var, time_integrator->getCurrentContext());
        const int Q_in_scr_idx =
            var_db->registerVariableAndContext(Q_in_var, var_db->getContext("SCRATCH"), IntVector<NDIM>(2));
        const int ls_in_idx = var_db->mapVariableAndContextToIndex(ls_in_var, time_integrator->getCurrentContext());
        const int ls_out_idx = var_db->mapVariableAndContextToIndex(ls_out_var, time_integrator->getCurrentContext());
        const int area_in_idx = var_db->mapVariableAndContextToIndex(time_integrator->getAreaVariable(ls_in_var),
                                                                     time_integrator->getCurrentContext());
        out_bdry_oper->registerAreaAndLSInsideIndex(area_in_idx, ls_in_idx);

        const int vol_in_idx = var_db->mapVariableAndContextToIndex(time_integrator->getVolumeVariable(ls_in_var),
                                                                    time_integrator->getCurrentContext());
        const int vol_out_idx = var_db->mapVariableAndContextToIndex(time_integrator->getVolumeVariable(ls_out_var),
                                                                     time_integrator->getCurrentContext());

        // Close the restart manager.
        RestartManager::getManager()->closeRestartFile();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        double dt = time_integrator->getMaximumTimeStepSize();

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            time_integrator->allocatePatchData(u_draw_idx, loop_time);
            u_fcn->setDataOnPatchHierarchy(u_draw_idx, u_draw_var, patch_hierarchy, loop_time);
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            time_integrator->deallocatePatchData(u_draw_idx);
            calculateTotalAmounts(patch_hierarchy, loop_time, Q_in_idx, Q_out_idx, vol_in_idx, vol_out_idx);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                time_integrator->allocatePatchData(u_draw_idx, loop_time);
                u_fcn->setDataOnPatchHierarchy(u_draw_idx, u_draw_var, patch_hierarchy, loop_time);
                if (output_bdry_info)
                {
                    outputBdryInfo(Q_in_idx,
                                   Q_in_scr_idx,
                                   ls_in_idx,
                                   vol_in_idx,
                                   area_in_idx,
                                   loop_time,
                                   iteration_num,
                                   patch_hierarchy,
                                   "in_bdry_info_");
                    outputBdryInfo(Q_out_idx,
                                   Q_out_scr_idx,
                                   ls_out_idx,
                                   vol_out_idx,
                                   area_in_idx,
                                   loop_time,
                                   iteration_num,
                                   patch_hierarchy,
                                   "out_bdry_info_");
                }
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                time_integrator->deallocatePatchData(u_draw_idx);
                calculateTotalAmounts(patch_hierarchy, loop_time, Q_in_idx, Q_out_idx, vol_in_idx, vol_out_idx);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % dump_postproc_interval == 0 || last_step))
            {
                postprocess_data(patch_hierarchy,
                                 time_integrator,
                                 Q_in_var,
                                 Q_out_var,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);
            }
        }

        if (!periodic_domain) delete Q_in_bcs[0];
        if (!periodic_domain) delete Q_out_bcs[0];
    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main

void
postprocess_data(Pointer<PatchHierarchy<NDIM>> hierarchy,
                 Pointer<LSAdvDiffIntegrator> integrator,
                 Pointer<CellVariable<NDIM, double>> Q_in_var,
                 Pointer<CellVariable<NDIM, double>> Q_out_var,
                 const int iteration_num,
                 const double loop_time,
                 const std::string& dirname)
{
    std::string file_name = dirname + "/hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    ComponentSelector hier_data;
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(Q_in_var, integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(Q_out_var, integrator->getCurrentContext()));
    hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
}

void
output_to_file(const int Q_idx,
               const int area_idx,
               const int vol_idx,
               const int ls_interp_idx,
               const std::string& file_name,
               const double loop_time,
               Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<hier::Variable<NDIM>> Q_var;
    var_db->mapIndexToVariable(Q_idx, Q_var);
    std::ofstream bdry_stream;
    if (SAMRAI_MPI::getRank() == 0) bdry_stream.open(file_name.c_str(), std::ofstream::out);
    // data structure to hold bdry data : (theta, bdry_val)
    std::vector<double> theta_data, val_data;
    // We only care about data on the finest level
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(hierarchy->getFinestLevelNumber());
    double integral = 0.0;
    double tot_area = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
        Pointer<CellData<NDIM, double>> area_data = patch->getPatchData(area_idx);
        Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
        Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(ls_interp_idx);
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        const double* const x_low = pgeom->getXLower();
        const Box<NDIM>& box = patch->getBox();
        const hier::Index<NDIM>& idx_low = box.lower();
        for (CellIterator<NDIM> ci(box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            const double area = (*area_data)(idx);
            if (area > 0.0 && (*vol_data)(idx) > 0.0)
            {
                std::array<VectorNd, 2> X_bounds;
                int l = 0;
                const NodeIndex<NDIM> idx_ll(idx, IntVector<NDIM>(0, 0)), idx_uu(idx, IntVector<NDIM>(1, 1)),
                    idx_lu(idx, IntVector<NDIM>(0, 1)), idx_ul(idx, IntVector<NDIM>(1, 0));
                const double phi_ll = (*ls_data)(idx_ll), phi_uu = (*ls_data)(idx_uu), phi_lu = (*ls_data)(idx_lu),
                             phi_ul = (*ls_data)(idx_ul);
                VectorNd X_ll(idx(0), idx(1)), X_uu(idx(0) + 1.0, idx(1) + 1.0), X_lu(idx(0), idx(1) + 1.0),
                    X_ul(idx(0) + 1.0, idx(1));
                if (phi_ll * phi_lu < 0.0)
                {
                    X_bounds[l] = midpoint_value(X_ll, phi_ll, X_lu, phi_lu);
                    l++;
                }
                if (phi_lu * phi_uu < 0.0)
                {
                    X_bounds[l] = midpoint_value(X_lu, phi_lu, X_uu, phi_uu);
                    l++;
                }
                if (phi_uu * phi_ul < 0.0)
                {
                    X_bounds[l] = midpoint_value(X_uu, phi_uu, X_ul, phi_ul);
                    l++;
                }
                if (phi_ul * phi_ll < 0.0)
                {
                    X_bounds[l] = midpoint_value(X_ul, phi_ul, X_ll, phi_ll);
                    l++;
                }
                TBOX_ASSERT(l == 2);
                VectorNd X = 0.5 * (X_bounds[0] + X_bounds[1]);
                VectorNd X_phys;
                for (int d = 0; d < NDIM; ++d) X_phys[d] = x_low[d] + dx[d] * (X(d) - static_cast<double>(idx_low(d)));
                double cur_theta = std::atan2(X_phys[1], X_phys[0]);
                double cur_r = X_phys.norm();
                double ref_theta =
                    cur_theta + (a * cur_r + b / cur_r) * (7.5 / (2.0 * M_PI)) * std::cos(2 * M_PI * loop_time / 7.5) -
                    (a * cur_r + b / cur_r) * (7.5 / (2.0 * M_PI));
                VectorNd X_new_coords = { cur_r * std::cos(ref_theta), cur_r * std::sin(ref_theta) };
                X_new_coords[0] -= 1.521;
                X_new_coords[1] -= 1.503;
                double theta = std::atan2(X_new_coords[1], X_new_coords[0]);
                theta_data.push_back(theta);
                // Calculate a delta theta for integral calculations
                double d_theta = area;

                // Do least squares linear approximation to find Q_val
                Box<NDIM> box_ls(idx, idx);
                box_ls.grow(2);
                std::vector<double> Q_vals;
                std::vector<VectorNd> X_vals;
                for (CellIterator<NDIM> ci(box_ls); ci; ci++)
                {
                    const CellIndex<NDIM>& idx_c = ci();
                    if ((*vol_data)(idx_c) > 0.0)
                    {
                        // Use this point
                        Q_vals.push_back((*Q_data)(idx_c));
                        X_vals.push_back(find_cell_centroid(idx_c, *ls_data));
                    }
                }
                const int m = Q_vals.size();
                MatrixXd A(MatrixXd::Zero(m, NDIM + 1));
                VectorXd U(VectorXd::Zero(m));
                for (size_t i = 0; i < Q_vals.size(); ++i)
                {
                    const VectorNd disp = X_vals[i] - X;
                    double w = std::sqrt(std::exp(-disp.norm() * disp.norm()));
                    A(i, 2) = w * disp[1];
                    A(i, 1) = w * disp[0];
                    A(i, 0) = w;
                    U(i) = w * Q_vals[i];
                }

                VectorXd soln = A.fullPivHouseholderQr().solve(U);
                val_data.push_back(soln(0));
                integral += soln(0) * d_theta;
                tot_area += d_theta;
            }
        }
    }
    integral = SAMRAI_MPI::sumReduction(integral);
    tot_area = SAMRAI_MPI::sumReduction(tot_area);
    pout << "Integral at time: " << loop_time << " for variable: " << Q_var->getName()
         << " is: " << std::setprecision(12) << integral << "\n";
    pout << "Area     at time: " << loop_time << " for variable: " << Q_var->getName()
         << " is: " << std::setprecision(12) << tot_area << "\n";
    // Now we need to send the data to processor rank 0 for outputting
    if (SAMRAI_MPI::getRank() == 0)
    {
        const int num_procs = SAMRAI_MPI::getNodes();
        std::vector<int> data_per_proc(num_procs - 1);
        for (int i = 1; i < num_procs; ++i)
        {
            MPI_Recv(&data_per_proc[i - 1], 1, MPI_INT, i, 0, SAMRAI_MPI::commWorld, nullptr);
        }
        std::vector<std::vector<double>> theta_per_proc(num_procs - 1), val_per_proc(num_procs - 1);
        for (int i = 1; i < num_procs; ++i)
        {
            theta_per_proc[i - 1].resize(data_per_proc[i - 1]);
            val_per_proc[i - 1].resize(data_per_proc[i - 1]);
            MPI_Recv(
                theta_per_proc[i - 1].data(), data_per_proc[i - 1], MPI_DOUBLE, i, 0, SAMRAI_MPI::commWorld, nullptr);
            MPI_Recv(
                val_per_proc[i - 1].data(), data_per_proc[i - 1], MPI_DOUBLE, i, 0, SAMRAI_MPI::commWorld, nullptr);
        }
        // Root processor now has all the data. Sort it and print it
        std::map<double, double> theta_val_data;
        // Start with root processor
        for (size_t i = 0; i < theta_data.size(); ++i) theta_val_data[theta_data[i]] = val_data[i];
        // Now loop through remaining processors
        for (int i = 1; i < num_procs; ++i)
        {
            for (size_t j = 0; j < theta_per_proc[i - 1].size(); ++j)
            {
                theta_val_data[theta_per_proc[i - 1][j]] = val_per_proc[i - 1][j];
            }
        }
        bdry_stream << std::setprecision(10) << loop_time << "\n";
        for (const auto& theta_val_pair : theta_val_data)
        {
            bdry_stream << theta_val_pair.first << " " << theta_val_pair.second << "\n";
        }
        bdry_stream.close();
    }
    else
    {
        TBOX_ASSERT(theta_data.size() == val_data.size());
        int num_data = theta_data.size();
        MPI_Send(&num_data, 1, MPI_INT, 0, 0, SAMRAI_MPI::commWorld);
        MPI_Send(theta_data.data(), theta_data.size(), MPI_DOUBLE, 0, 0, SAMRAI_MPI::commWorld);
        MPI_Send(val_data.data(), val_data.size(), MPI_DOUBLE, 0, 0, SAMRAI_MPI::commWorld);
    }
}
