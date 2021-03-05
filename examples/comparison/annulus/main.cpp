#include <IBAMR_config.h>
#include <IBTK_config.h>

#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/app_namespaces.h>

#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/muParserCartGridFunction.h"
#include "ibtk/muParserRobinBcCoefs.h"
#include <ibtk/AppInitializer.h>

#include "LS/LSCutCellLaplaceOperator.h"
#include "LS/LSFromLevelSet.h"
#include "LS/QInitial.h"
#include "LS/SemiLagrangianAdvIntegrator.h"

#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "LoadBalancer.h"
#include "SAMRAI_config.h"
#include "StandardTagAndInitialize.h"
#include "tbox/Pointer.h"

#include <petscsys.h>

#include <Eigen/Dense>

#include <memory>
#include <utility>

// Local includes
#include "LSFcn.h"
#include "QFcn.h"
#include "RadialBoundaryCond.h"

using namespace LS;

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
        Pointer<SemiLagrangianAdvIntegrator> time_integrator = new SemiLagrangianAdvIntegrator(
            "SemiLagrangianAdvIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"),
            false);

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

        // Setup advected quantity
        Pointer<CellVariable<NDIM, double>> Q_var = new CellVariable<NDIM, double>("Q");
        Pointer<QFcn> Q_init = new QFcn("QInit", grid_geometry, app_initializer->getComponentDatabase("QInitial"));
        const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;
        std::vector<RobinBcCoefStrategy<NDIM>*> Q_bcs(1);
        if (!periodic_domain)
        {
            const std::string Q_bcs_name = "Q_bcs";
            Q_bcs[0] =
                new muParserRobinBcCoefs(Q_bcs_name, app_initializer->getComponentDatabase(Q_bcs_name), grid_geometry);
        }

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
        Pointer<NodeVariable<NDIM, double>> ls_var = new NodeVariable<NDIM, double>("LS");
        time_integrator->registerLevelSetVariable(ls_var);
        bool use_ls_fcn = input_db->getBool("USING_LS_FCN");
        Pointer<LSFcn> ls_fcn = new LSFcn("LSFcn", app_initializer->getComponentDatabase("LSFcn"));
        Pointer<LSFromLevelSet> vol_fcn = new LSFromLevelSet("LSFromLevelSet", patch_hierarchy);
        vol_fcn->registerLSFcn(ls_fcn);
        time_integrator->registerLevelSetVolFunction(ls_var, vol_fcn);

        // Setup forcing function
        Pointer<CellVariable<NDIM, double>> F_var = new CellVariable<NDIM, double>("F");
        Pointer<CartGridFunction> F_fcn =
            new muParserCartGridFunction("F_Fcn", app_initializer->getComponentDatabase("F_Fcn"), grid_geometry);
        time_integrator->registerSourceTerm(F_var);
        time_integrator->setSourceTermFunction(F_var, F_fcn);

        time_integrator->registerTransportedQuantity(Q_var);
        time_integrator->setAdvectionVelocity(Q_var, u_var);
        time_integrator->setInitialConditions(Q_var, Q_init);
        time_integrator->setPhysicalBcCoef(Q_var, Q_bcs[0]);
        time_integrator->setDiffusionCoefficient(Q_var, input_db->getDoubleWithDefault("D_coef", 0.0));
        time_integrator->restrictToLevelSet(Q_var, ls_var);
        time_integrator->setSourceTerm(Q_var, F_var);

        // Set up diffusion operators
        Pointer<RadialBoundaryCond> rhs_bdry_oper =
            new RadialBoundaryCond("RHSBdryOperator", app_initializer->getComponentDatabase("BdryOperator"));
        Pointer<RadialBoundaryCond> sol_bdry_oper =
            new RadialBoundaryCond("SOLBdryOperator", app_initializer->getComponentDatabase("BdryOperator"));
        Pointer<LSCutCellLaplaceOperator> rhs_oper = new LSCutCellLaplaceOperator(
            "LSCutCellRHSOperator", app_initializer->getComponentDatabase("LSCutCellOperator"), false);
        Pointer<LSCutCellLaplaceOperator> sol_oper = new LSCutCellLaplaceOperator(
            "LSCutCellOperator", app_initializer->getComponentDatabase("LSCutCellOperator"), false);
        rhs_oper->setBoundaryConditionOperator(rhs_bdry_oper);
        sol_oper->setBoundaryConditionOperator(sol_bdry_oper);
        time_integrator->setHelmholtzRHSOperator(Q_var, rhs_oper);
        Pointer<PETScKrylovPoissonSolver> Q_helmholtz_solver = new PETScKrylovPoissonSolver(
            "PoissonSolver", app_initializer->getComponentDatabase("PoissonSolver"), "poisson_solve_");
        Q_helmholtz_solver->setOperator(sol_oper);
        Q_helmholtz_solver->setNullspace(true);
        time_integrator->setHelmholtzSolver(Q_var, Q_helmholtz_solver);

        // Create a scratch volume, area, and level set
        bool output_bdry_info = input_db->getBool("OUTPUT_BDRY_INFO");
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int ls_idx =
            var_db->registerVariableAndContext(ls_var, var_db->getContext("SCRATCH"), IntVector<NDIM>(4));
        Pointer<CellVariable<NDIM, double>> vol_var = new CellVariable<NDIM, double>("vol");
        Pointer<CellVariable<NDIM, double>> area_var = new CellVariable<NDIM, double>("area");
        const int vol_idx =
            var_db->registerVariableAndContext(vol_var, var_db->getContext("SCRATCH"), IntVector<NDIM>(4));
        const int area_idx =
            var_db->registerVariableAndContext(area_var, var_db->getContext("SCRATCH"), IntVector<NDIM>(4));

        // Set up visualization plot file writer.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Register a drawing variable with the data writer
        Pointer<CellVariable<NDIM, double>> u_draw_var = new CellVariable<NDIM, double>("U", NDIM);
        const int u_draw_idx = var_db->registerVariableAndContext(u_draw_var, var_db->getContext("Draw"));
        visit_data_writer->registerPlotQuantity("U", "VECTOR", u_draw_idx);
        const int Q_exact_idx = var_db->registerVariableAndContext(Q_var, var_db->getContext("Draw"));
        bool draw_exact = input_db->getBool("DRAW_EXACT");
        if (draw_exact) visit_data_writer->registerPlotQuantity("Exact", "SCALAR", Q_exact_idx);
        if (draw_exact) visit_data_writer->registerPlotQuantity("Exact LS", "SCALAR", ls_idx);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, time_integrator->getCurrentContext());
        const int Q_scr_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("SCRATCH"), IntVector<NDIM>(1));

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
            if (draw_exact) time_integrator->allocatePatchData(Q_exact_idx, loop_time);
            if (draw_exact) Q_init->setDataOnPatchHierarchy(Q_exact_idx, Q_var, patch_hierarchy, loop_time);
            u_fcn->setDataOnPatchHierarchy(u_draw_idx, u_draw_var, patch_hierarchy, loop_time);
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            time_integrator->deallocatePatchData(u_draw_idx);
            if (draw_exact) time_integrator->deallocatePatchData(Q_exact_idx);
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
                if (draw_exact)
                {
                    time_integrator->allocatePatchData(Q_exact_idx, loop_time);
                    time_integrator->allocatePatchData(ls_idx, loop_time);
                    time_integrator->allocatePatchData(vol_idx, loop_time);
                    vol_fcn->updateVolumeAreaSideLS(
                        vol_idx, vol_var, -1, nullptr, -1, nullptr, ls_idx, ls_var, loop_time, false);
                    Q_init->setLSIndex(ls_idx, vol_idx);
                    Q_init->setDataOnPatchHierarchy(Q_exact_idx, Q_var, patch_hierarchy, loop_time);
                }
                u_fcn->setDataOnPatchHierarchy(u_draw_idx, u_draw_var, patch_hierarchy, loop_time);
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                time_integrator->deallocatePatchData(u_draw_idx);
                if (draw_exact)
                {
                    time_integrator->deallocatePatchData(Q_exact_idx);
                    time_integrator->deallocatePatchData(ls_idx);
                    time_integrator->deallocatePatchData(vol_idx);
                }
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
        }

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        const int Q_err_idx = var_db->registerVariableAndContext(Q_var, var_db->getContext("Error Context"));

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_err_idx, loop_time);
            level->allocatePatchData(ls_idx, loop_time);
            level->allocatePatchData(vol_idx, loop_time);
        }

        ls_fcn->setDataOnPatchHierarchy(ls_idx, ls_var, patch_hierarchy, loop_time);
        vol_fcn->updateVolumeAreaSideLS(vol_idx, vol_var, -1, nullptr, -1, nullptr, ls_idx, ls_var, loop_time, false);
        Q_init->setLSIndex(ls_idx, vol_idx);

        Q_init->setDataOnPatchHierarchy(Q_err_idx, Q_var, patch_hierarchy, loop_time);

        Pointer<HierarchyMathOps> hier_math_ops = new HierarchyMathOps("HierarchyMathOps", patch_hierarchy);
        const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());

                Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);
                Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    (*wgt_data)(idx) *= (*vol_data)(idx);
                }
            }
        }

        pout << "Norms of exact solution at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(Q_err_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(Q_err_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(Q_err_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        hier_cc_data_ops.subtract(Q_idx, Q_idx, Q_err_idx);
        pout << "Error in " << Q_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(Q_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(Q_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(Q_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);
                Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    (*wgt_data)(idx) *= (*vol_data)(idx) < 1.0 ? 0.0 : 1.0;
                }
            }
        }

        pout << "Error without cut cells in " << Q_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(Q_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(Q_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(Q_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        if (dump_viz_data)
        {
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                time_integrator->allocatePatchData(u_draw_idx, loop_time);
                if (draw_exact) time_integrator->allocatePatchData(Q_exact_idx, loop_time);
                if (draw_exact) Q_init->setDataOnPatchHierarchy(Q_exact_idx, Q_var, patch_hierarchy, loop_time);
                u_fcn->setDataOnPatchHierarchy(u_draw_idx, u_draw_var, patch_hierarchy, loop_time);
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
                time_integrator->deallocatePatchData(u_draw_idx);
                if (draw_exact) time_integrator->deallocatePatchData(Q_exact_idx);
            }
        }

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_err_idx);
        }

        if (!periodic_domain) delete Q_bcs[0];
    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main
