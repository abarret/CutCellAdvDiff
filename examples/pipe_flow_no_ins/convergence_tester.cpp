// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// GENERAL CONFIGURATION

#include <SAMRAI_config.h>

#include <petsc.h>

#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>

#include <libmesh/equation_systems.h>
#include <libmesh/exact_solution.h>
#include <libmesh/mesh.h>
#include <libmesh/exodusII_io.h>

#include <ibamr/app_namespaces.h>

#include "ibtk/AppInitializer.h"
#include "ibtk/IBTKInit.h"

#include "LS/SemiLagrangianAdvIntegrator.h"

#include <tbox/Database.h>
#include <tbox/HDFDatabase.h>
#include <tbox/InputDatabase.h>
#include <tbox/InputManager.h>
#include <tbox/MathUtilities.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/SAMRAIManager.h>
#include <tbox/Utilities.h>

#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <ComponentSelector.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <PatchHierarchy.h>
#include <SideVariable.h>
#include <VariableDatabase.h>
#include <VisItDataWriter.h>

#include "ibamr/app_namespaces.h"

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    {
        if (argc != 2)
        {
            tbox::pout << "USAGE:  " << argv[0] << " <input filename>\n"
                       << "  options:\n"
                       << "  none at this time" << endl;
            TBOX_ERROR("Not enough input arguments.\n");
            return (-1);
        }

        string input_filename = argv[1];
        tbox::plog << "input_filename = " << input_filename << endl;

        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        int coarse_hier_dump_interval = 0;
        int fine_hier_dump_interval = 0;
        if (main_db->keyExists("hier_dump_interval"))
        {
            coarse_hier_dump_interval = main_db->getInteger("hier_dump_interval");
            fine_hier_dump_interval = main_db->getInteger("hier_dump_interval");
        }
        else if (main_db->keyExists("coarse_hier_dump_interval") && main_db->keyExists("fine_hier_dump_interval"))
        {
            coarse_hier_dump_interval = main_db->getInteger("coarse_hier_dump_interval");
            fine_hier_dump_interval = main_db->getInteger("fine_hier_dump_interval");
        }
        else
        {
            TBOX_ERROR("hierarchy dump intervals not specified in input file. . .\n");
        }

        string coarse_hier_dump_dirname;
        if (main_db->keyExists("coarse_hier_dump_dirname"))
        {
            coarse_hier_dump_dirname = main_db->getString("coarse_hier_dump_dirname");
        }
        else
        {
            TBOX_ERROR("key `coarse_hier_dump_dirname' not specified in input file");
        }

        string fine_hier_dump_dirname;
        if (main_db->keyExists("fine_hier_dump_dirname"))
        {
            fine_hier_dump_dirname = main_db->getString("fine_hier_dump_dirname");
        }
        else
        {
            TBOX_ERROR("key `fine_hier_dump_dirname' not specified in input file");
        }

        // Create major algorithm and data objects which comprise application.
        tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geom =
            new geom::CartesianGridGeometry<NDIM>("CartesianGeometry", input_db->getDatabase("CartesianGeometry"));

        // Initialize variables.
        hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();

        tbox::Pointer<hier::VariableContext> current_ctx = var_db->getContext("SemiLagrangianAdvIntegrator::CURRENT");
        tbox::Pointer<hier::VariableContext> scratch_ctx = var_db->getContext("SemiLagrangianAdvIntegrator::SCRATCH");

        tbox::Pointer<pdat::CellVariable<NDIM, double>> Q_var = new pdat::CellVariable<NDIM, double>("Q_in");
        const int Q_idx = var_db->registerVariableAndContext(Q_var, current_ctx);
        const int Q_interp_idx = var_db->registerClonedPatchDataIndex(Q_var, Q_idx);
        const int Q_scr_idx = var_db->registerVariableAndContext(Q_var, scratch_ctx, 2);

        tbox::Pointer<pdat::CellVariable<NDIM, double>> vol_var = new pdat::CellVariable<NDIM, double>("LS_VolVar");
        const int vol_idx = var_db->registerVariableAndContext(vol_var, current_ctx);

        tbox::Pointer<pdat::NodeVariable<NDIM, double>> ls_var = new pdat::NodeVariable<NDIM, double>("LS");
        const int ls_idx = var_db->registerVariableAndContext(ls_var, current_ctx);

        // Set up visualization plot file writer.
        tbox::Pointer<appu::VisItDataWriter<NDIM>> visit_data_writer =
            new appu::VisItDataWriter<NDIM>("VisIt Writer", main_db->getString("viz_dump_dirname"), 1);
        visit_data_writer->registerPlotQuantity("Q", "SCALAR", Q_idx);
        visit_data_writer->registerPlotQuantity("Q interp", "SCALAR", Q_interp_idx);
        visit_data_writer->registerPlotQuantity("Volume", "SCALAR", vol_idx);
        visit_data_writer->registerPlotQuantity("LS", "SCALAR", ls_idx);

        // Time step loop.
        double loop_time = 0.0;
        int coarse_iteration_num = coarse_hier_dump_interval;
        int fine_iteration_num = fine_hier_dump_interval;

        bool files_exist = true;
        for (; files_exist;
             coarse_iteration_num += coarse_hier_dump_interval, fine_iteration_num += fine_hier_dump_interval)
        {
            char temp_buf[128];

            sprintf(temp_buf, "%05d.samrai.%05d", coarse_iteration_num, IBTK_MPI::getRank());
            string coarse_file_name = coarse_hier_dump_dirname + "/" + "hier_data.";
            coarse_file_name += temp_buf;

            sprintf(temp_buf, "%05d.samrai.%05d", fine_iteration_num, IBTK_MPI::getRank());
            string fine_file_name = fine_hier_dump_dirname + "/" + "hier_data.";
            fine_file_name += temp_buf;

            for (int rank = 0; rank < IBTK_MPI::getNodes(); ++rank)
            {
                if (rank == IBTK_MPI::getRank())
                {
                    fstream coarse_fin, fine_fin;
                    coarse_fin.open(coarse_file_name.c_str(), ios::in);
                    fine_fin.open(fine_file_name.c_str(), ios::in);
                    if (!coarse_fin.is_open() || !fine_fin.is_open())
                    {
                        files_exist = false;
                    }
                    coarse_fin.close();
                    fine_fin.close();
                }
                IBTK_MPI::barrier();
            }

            if (!files_exist) break;

            tbox::pout << endl;
            tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << "processing data" << endl;
            tbox::pout << "     coarse iteration number = " << coarse_iteration_num << endl;
            tbox::pout << "     fine iteration number = " << fine_iteration_num << endl;
            tbox::pout << "     coarse file name = " << coarse_file_name << endl;
            tbox::pout << "     fine file name = " << fine_file_name << endl;

            // Read in data to post-process.
            hier::ComponentSelector hier_data;
            hier_data.setFlag(Q_idx);
            hier_data.setFlag(vol_idx);
            hier_data.setFlag(ls_idx);

            tbox::Pointer<tbox::HDFDatabase> coarse_hier_db = new tbox::HDFDatabase("coarse_hier_db");
            coarse_hier_db->open(coarse_file_name);

            tbox::Pointer<hier::PatchHierarchy<NDIM>> coarse_patch_hierarchy =
                new hier::PatchHierarchy<NDIM>("CoarsePatchHierarchy", grid_geom, false);
            coarse_patch_hierarchy->getFromDatabase(coarse_hier_db->getDatabase("PatchHierarchy"), hier_data);

            const double coarse_loop_time = coarse_hier_db->getDouble("loop_time");

            coarse_hier_db->close();

            tbox::Pointer<tbox::HDFDatabase> fine_hier_db = new tbox::HDFDatabase("fine_hier_db");
            fine_hier_db->open(fine_file_name);

            tbox::Pointer<hier::PatchHierarchy<NDIM>> fine_patch_hierarchy = new hier::PatchHierarchy<NDIM>(
                "FinePatchHierarchy", grid_geom->makeRefinedGridGeometry("FineGridGeometry", 2, false), false);
            fine_patch_hierarchy->getFromDatabase(fine_hier_db->getDatabase("PatchHierarchy"), hier_data);

            const double fine_loop_time = fine_hier_db->getDouble("loop_time");

            fine_hier_db->close();

            TBOX_ASSERT(tbox::MathUtilities<double>::equalEps(coarse_loop_time, fine_loop_time));
            loop_time = fine_loop_time;
            tbox::pout << "     loop time = " << loop_time << endl;

            tbox::Pointer<hier::PatchHierarchy<NDIM>> coarsened_fine_patch_hierarchy =
                fine_patch_hierarchy->makeCoarsenedPatchHierarchy("CoarsenedFinePatchHierarchy", 2, false);

            // Setup hierarchy operations objects.
            math::HierarchyCellDataOpsReal<NDIM, double> coarse_hier_cc_data_ops(
                coarse_patch_hierarchy, 0, coarse_patch_hierarchy->getFinestLevelNumber());
            HierarchyMathOps hier_math_ops("hier_math_ops", coarse_patch_hierarchy);
            hier_math_ops.setPatchHierarchy(coarse_patch_hierarchy);
            hier_math_ops.resetLevels(0, coarse_patch_hierarchy->getFinestLevelNumber());
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

            const IntVector<NDIM> cut_cell_width(input_db->getInteger("SKIP_CUT_CELLS"));
            for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = coarse_patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);
                    Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
                    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                    const double* const dx = pgeom->getDx();
                    for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = ci();
                        Box<NDIM> box(idx, idx);
                        box.grow(cut_cell_width);
                        bool next_to_cut = false;
                        for (CellIterator<NDIM> cii(box); cii; cii++)
                        {
                            const CellIndex<NDIM>& idx_i = cii();
                            if ((*vol_data)(idx_i) < 1.0) next_to_cut = true;
                        }
                        (*wgt_data)(idx) *= next_to_cut ? 0.0 : (*vol_data)(idx);
                    }
                }
            }

            // Allocate patch data.
            for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> level = coarse_patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(Q_interp_idx, loop_time);
                level->allocatePatchData(Q_scr_idx, loop_time);
            }

            for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> level = fine_patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(Q_interp_idx, loop_time);
                level->allocatePatchData(Q_scr_idx, loop_time);
            }

            for (int ln = 0; ln <= coarsened_fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(Q_idx, loop_time);
                level->allocatePatchData(Q_interp_idx, loop_time);
                level->allocatePatchData(Q_scr_idx, loop_time);
            }

            // Synchronize the coarse hierarchy data.
            for (int ln = coarse_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> coarser_level = coarse_patch_hierarchy->getPatchLevel(ln - 1);
                tbox::Pointer<hier::PatchLevel<NDIM>> finer_level = coarse_patch_hierarchy->getPatchLevel(ln);

                xfer::CoarsenAlgorithm<NDIM> coarsen_alg;
                tbox::Pointer<xfer::CoarsenOperator<NDIM>> coarsen_op;

                coarsen_op = grid_geom->lookupCoarsenOperator(Q_var, "CONSERVATIVE_COARSEN");
                coarsen_alg.registerCoarsen(Q_idx, Q_idx, coarsen_op);

                coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
            }

            // Synchronize the fine hierarchy data.
            for (int ln = fine_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> coarser_level = fine_patch_hierarchy->getPatchLevel(ln - 1);
                tbox::Pointer<hier::PatchLevel<NDIM>> finer_level = fine_patch_hierarchy->getPatchLevel(ln);

                xfer::CoarsenAlgorithm<NDIM> coarsen_alg;
                tbox::Pointer<xfer::CoarsenOperator<NDIM>> coarsen_op;

                coarsen_op = grid_geom->lookupCoarsenOperator(Q_var, "CONSERVATIVE_COARSEN");
                coarsen_alg.registerCoarsen(Q_idx, Q_idx, coarsen_op);

                coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
            }

            // Coarsen data from the fine hierarchy to the coarsened fine hierarchy.
            for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> dst_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
                tbox::Pointer<hier::PatchLevel<NDIM>> src_level = fine_patch_hierarchy->getPatchLevel(ln);

                tbox::Pointer<xfer::CoarsenOperator<NDIM>> coarsen_op;
                for (hier::PatchLevel<NDIM>::Iterator p(dst_level); p; p++)
                {
                    tbox::Pointer<hier::Patch<NDIM>> dst_patch = dst_level->getPatch(p());
                    tbox::Pointer<hier::Patch<NDIM>> src_patch = src_level->getPatch(p());
                    const hier::Box<NDIM>& coarse_box = dst_patch->getBox();
                    TBOX_ASSERT(hier::Box<NDIM>::coarsen(src_patch->getBox(), 2) == coarse_box);

                    coarsen_op = grid_geom->lookupCoarsenOperator(Q_var, "CONSERVATIVE_COARSEN");
                    coarsen_op->coarsen(*dst_patch, *src_patch, Q_interp_idx, Q_idx, coarse_box, 2);
                }
            }

            // Interpolate and copy data from the coarsened fine patch hierarchy to
            // the coarse patch hierarchy.
            for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> dst_level = coarse_patch_hierarchy->getPatchLevel(ln);
                tbox::Pointer<hier::PatchLevel<NDIM>> src_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);

                xfer::RefineAlgorithm<NDIM> refine_alg;
                tbox::Pointer<xfer::RefineOperator<NDIM>> refine_op;

                refine_op = grid_geom->lookupRefineOperator(Q_var, "CONSERVATIVE_LINEAR_REFINE");
                refine_alg.registerRefine(Q_interp_idx, Q_interp_idx, Q_scr_idx, refine_op);

                hier::ComponentSelector data_indices;
                data_indices.setFlag(Q_scr_idx);
                CartExtrapPhysBdryOp bc_helper(data_indices, "LINEAR");

                refine_alg.createSchedule(dst_level, src_level, ln - 1, coarse_patch_hierarchy, &bc_helper)
                    ->fillData(loop_time);
            }

            // Output plot data before taking norms of differences.
            visit_data_writer->writePlotData(coarse_patch_hierarchy, coarse_iteration_num, loop_time);

            // Compute norms of differences.
            coarse_hier_cc_data_ops.subtract(Q_interp_idx, Q_idx, Q_interp_idx);

            tbox::pout << "\n"
                       << "Error in " << Q_var->getName() << " at time " << loop_time << ":\n"
                       << "  Q_in_L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(Q_interp_idx, wgt_cc_idx) << "\n"
                       << "  Q_in_L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(Q_interp_idx, wgt_cc_idx) << "\n"
                       << "  Q_in_max-norm: " << coarse_hier_cc_data_ops.maxNorm(Q_interp_idx, wgt_cc_idx) << "\n";

            // Output plot data after taking norms of differences.
            visit_data_writer->writePlotData(coarse_patch_hierarchy, coarse_iteration_num + 1, loop_time);

            // Do the same thing for the FE data.
            string file_name;

            Mesh mesh_coarse(init.comm(), NDIM);
            file_name = coarse_hier_dump_dirname + "/" + "fe_mesh.";
            sprintf(temp_buf, "%05d", coarse_iteration_num);
            file_name += temp_buf;
            file_name += ".xda";
            mesh_coarse.read(file_name);

            Mesh mesh_fine(init.comm(), NDIM);
            file_name = fine_hier_dump_dirname + "/" + "fe_mesh.";
            sprintf(temp_buf, "%05d", fine_iteration_num);
            file_name += temp_buf;
            file_name += ".xda";
            mesh_fine.read(file_name);

            EquationSystems equation_systems_coarse(mesh_coarse);
            file_name = coarse_hier_dump_dirname + "/" + "fe_equation_systems.";
            sprintf(temp_buf, "%05d", coarse_iteration_num);
            file_name += temp_buf;
            equation_systems_coarse.read(
                file_name,
                (EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA));

            EquationSystems equation_systems_fine(mesh_fine);
            file_name = fine_hier_dump_dirname + "/" + "fe_equation_systems.";
            sprintf(temp_buf, "%05d", fine_iteration_num);
            file_name += temp_buf;
            equation_systems_fine.read(
                file_name,
                (EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA));

            ExactSolution error_estimator(equation_systems_coarse);
            error_estimator.attach_reference_solution(&equation_systems_fine);

            error_estimator.compute_error("SurfaceConcentration", "SurfaceConcentration");
            double sf_error[3];
            sf_error[0] = error_estimator.l1_error("SurfaceConcentration", "SurfaceConcentration");
            sf_error[1] = error_estimator.l2_error("SurfaceConcentration", "SurfaceConcentration");
            sf_error[2] = error_estimator.l_inf_error("SurfaceConcentration", "SurfaceConcentration");

            error_estimator.compute_error("Q_in", "Q_in");
            double fl_error[3];
            fl_error[0] = error_estimator.l1_error("Q_in", "Q_in");
            fl_error[1] = error_estimator.l2_error("Q_in", "Q_in");
            fl_error[2] = error_estimator.l_inf_error("Q_in", "Q_in");
            tbox::pout << "\n"
                       << "Error in Surface Concentration at time " << loop_time << ":\n"
                       << "  Surface Concentration_L1-norm:  " << sf_error[0] << "\n"
                       << "  Surface Concentration_L2-norm:  " << sf_error[1] << "\n"
                       << "  Surface Concentration_max-norm: " << sf_error[2] << "\n";
            tbox::pout << "\n"
                       << "Error in Q on surface at time " << loop_time << ":\n"
                       << "  Q on surface_L1-norm:  " << fl_error[0] << "\n"
                       << "  Q on surface_L2-norm:  " << fl_error[1] << "\n"
                       << "  Q on surface_max-norm: " << fl_error[2] << "\n";

            tbox::pout << endl;
            tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << endl;
        }
    }
    return 0;
} // main
