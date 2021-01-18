// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
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
#include <ibamr/config.h>

#include <petscsys.h>

#include <SAMRAI_config.h>

// IBAMR INCLUDES
#include <ibamr/app_namespaces.h>

// IBTK INCLUDES
#include <ibtk/AppInitializer.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

#include "LS/LSFindCellVolume.h"

#include "BoxArray.h"
#include "CartesianPatchGeometry.h"
#include "CoarseFineBoundary.h"
#include "LSFcn.h"
#include "PatchGeometry.h"
#include "RefineAlgorithm.h"

using namespace LS;

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name> [PETSc options]                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // Parse command line options, set some standard options from the input
    // file, and enable file logging.
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Retrieve "Main" section of the input database.
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

    // Create major algorithm and data objects that comprise application.
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    // Initialize variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    Pointer<VariableContext> current_ctx = var_db->getContext("SemiLagrangianAdvIntegrator::CURRENT");
    Pointer<VariableContext> scratch_ctx = var_db->getContext("SemiLagrangianAdvIntegrator::SCRATCH");

    Pointer<CellVariable<NDIM, double>> Q_var = new CellVariable<NDIM, double>("Q");
    const int Q_idx = var_db->registerVariableAndContext(Q_var, current_ctx);
    const int Q_interp_idx = var_db->registerClonedPatchDataIndex(Q_var, Q_idx);
    const int Q_scr_idx = var_db->registerVariableAndContext(Q_var, scratch_ctx, 2);

    Pointer<NodeVariable<NDIM, double>> ls_var = new NodeVariable<NDIM, double>("LS_var");
    const int ls_idx = var_db->registerVariableAndContext(ls_var, current_ctx, IntVector<NDIM>(2));

    Pointer<CellVariable<NDIM, double>> vol_var = new CellVariable<NDIM, double>("Volume In");
    const int vol_idx = var_db->registerVariableAndContext(vol_var, current_ctx, IntVector<NDIM>(2));

    Pointer<CellVariable<NDIM, int>> rank_var = new CellVariable<NDIM, int>("Rank");
    const int rank_idx = var_db->registerVariableAndContext(rank_var, current_ctx);

    // Set up visualization plot file writer.
    Pointer<VisItDataWriter<NDIM>> visit_data_writer =
        new VisItDataWriter<NDIM>("VisIt Writer", main_db->getString("viz_dump_dirname"), 1);
    visit_data_writer->registerPlotQuantity("Q", "SCALAR", Q_idx);
    visit_data_writer->registerPlotQuantity("Q interp", "SCALAR", Q_interp_idx);
    visit_data_writer->registerPlotQuantity("ls", "SCALAR", ls_idx);
    visit_data_writer->registerPlotQuantity("vol", "SCALAR", vol_idx);
    visit_data_writer->registerPlotQuantity("rank", "SCALAR", rank_idx);

    // Time step loop.
    double loop_time = 0.0;
    int coarse_iteration_num = coarse_hier_dump_interval;
    int fine_iteration_num = fine_hier_dump_interval;

    bool files_exist = true;
    for (; files_exist;
         coarse_iteration_num += coarse_hier_dump_interval, fine_iteration_num += fine_hier_dump_interval)
    {
        char temp_buf[128];

        sprintf(temp_buf, "%05d.samrai.%05d", coarse_iteration_num, SAMRAI_MPI::getRank());
        string coarse_file_name = coarse_hier_dump_dirname + "/" + "hier_data.";
        coarse_file_name += temp_buf;

        sprintf(temp_buf, "%05d.samrai.%05d", fine_iteration_num, SAMRAI_MPI::getRank());
        string fine_file_name = fine_hier_dump_dirname + "/" + "hier_data.";
        fine_file_name += temp_buf;

        for (int rank = 0; rank < SAMRAI_MPI::getNodes(); ++rank)
        {
            if (rank == SAMRAI_MPI::getRank())
            {
                fstream coarse_fin, fine_fin;
                coarse_fin.open(coarse_file_name.c_str(), ios::in);
                fine_fin.open(fine_file_name.c_str(), ios::in);
                if (!coarse_fin.is_open() || !fine_fin.is_open())
                {
                    std::cout << "couldn't find file " << rank << "\n";
                    std::cout << "File name was: " << coarse_file_name << "\n and " << fine_file_name << "\n";
                    files_exist = false;
                }
                else
                {
                    std::cout << "found file " << rank << "\n";
                }
                coarse_fin.close();
                fine_fin.close();
            }
            SAMRAI_MPI::barrier();
        }

        if (!files_exist) break;

        pout << endl;
        pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout << "processing data" << endl;
        pout << "     coarse iteration number = " << coarse_iteration_num << endl;
        pout << "     fine iteration number = " << fine_iteration_num << endl;
        pout << "     coarse file name = " << coarse_file_name << endl;
        pout << "     fine file name = " << fine_file_name << endl;

        // Read in data to post-process.
        ComponentSelector hier_data;
        hier_data.setFlag(Q_idx);

        Pointer<HDFDatabase> coarse_hier_db = new HDFDatabase("coarse_hier_db");
        coarse_hier_db->open(coarse_file_name);

        Pointer<PatchHierarchy<NDIM>> coarse_patch_hierarchy =
            new PatchHierarchy<NDIM>("CoarsePatchHierarchy", grid_geom, false);
        coarse_patch_hierarchy->getFromDatabase(coarse_hier_db->getDatabase("PatchHierarchy"), hier_data);

        const double coarse_loop_time = coarse_hier_db->getDouble("loop_time");

        coarse_hier_db->close();

        Pointer<HDFDatabase> fine_hier_db = new HDFDatabase("fine_hier_db");
        fine_hier_db->open(fine_file_name);

        Pointer<PatchHierarchy<NDIM>> fine_patch_hierarchy = new PatchHierarchy<NDIM>(
            "FinePatchHierarchy", grid_geom->makeRefinedGridGeometry("FineGridGeometry", 2, false), false);
        fine_patch_hierarchy->getFromDatabase(fine_hier_db->getDatabase("PatchHierarchy"), hier_data);

        const double fine_loop_time = fine_hier_db->getDouble("loop_time");

        fine_hier_db->close();

        TBOX_ASSERT(MathUtilities<double>::equalEps(coarse_loop_time, fine_loop_time));
        loop_time = fine_loop_time;
        pout << "     loop time = " << loop_time << endl;

        Pointer<PatchHierarchy<NDIM>> coarsened_fine_patch_hierarchy =
            fine_patch_hierarchy->makeCoarsenedPatchHierarchy("CoarsenedFinePatchHierarchy", 2, false);

        // Setup hierarchy operations objects.
        HierarchyCellDataOpsReal<NDIM, double> coarse_hier_cc_data_ops(
            coarse_patch_hierarchy, 0, coarse_patch_hierarchy->getFinestLevelNumber());
        HierarchyMathOps hier_math_ops("hier_math_ops", coarse_patch_hierarchy);
        hier_math_ops.setPatchHierarchy(coarse_patch_hierarchy);
        hier_math_ops.resetLevels(0, coarse_patch_hierarchy->getFinestLevelNumber());
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        // Allocate patch data.
        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = coarse_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_interp_idx, loop_time);
            level->allocatePatchData(ls_idx, loop_time);
            level->allocatePatchData(vol_idx, loop_time);
            level->allocatePatchData(Q_scr_idx, loop_time);
            level->allocatePatchData(rank_idx, loop_time);
        }

        for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = fine_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_interp_idx, loop_time);
            level->allocatePatchData(Q_scr_idx, loop_time);
            level->allocatePatchData(ls_idx, loop_time);
            level->allocatePatchData(vol_idx, loop_time);
            level->allocatePatchData(rank_idx, loop_time);
        }

        for (int ln = 0; ln <= coarsened_fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_interp_idx, loop_time);
            level->allocatePatchData(Q_idx, loop_time);
            level->allocatePatchData(Q_scr_idx, loop_time);
            level->allocatePatchData(ls_idx, loop_time);
            level->allocatePatchData(vol_idx, loop_time);
            level->allocatePatchData(rank_idx, loop_time);
        }

        Pointer<CartGridFunction> ls_fcn = new LSFcn("SetLSValue", app_initializer->getComponentDatabase("LSFcn"));
        {
            // Coarse patch hierarchy
            Pointer<LSFindCellVolume> vol_fcn = new LSFindCellVolume("VolFcn", coarse_patch_hierarchy);
            ls_fcn->setDataOnPatchHierarchy(ls_idx, ls_var, coarse_patch_hierarchy, loop_time);

            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<ITC> ghost_cell_comps(1);
            ghost_cell_comps[0] = ITC(ls_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            HierarchyGhostCellInterpolation hier_ghost_cell;
            hier_ghost_cell.initializeOperatorState(ghost_cell_comps, coarse_patch_hierarchy);
            hier_ghost_cell.fillData(loop_time);

            vol_fcn->updateVolumeAndArea(vol_idx, vol_var, IBTK::invalid_index, nullptr, ls_idx, ls_var, true);
            for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = coarse_patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CellData<NDIM, int>> rank_data = patch->getPatchData(rank_idx);
                    rank_data->fillAll(SAMRAI_MPI::getRank());
                }
            }
        }

        {
            // Fine patch hierarchy
            Pointer<LSFindCellVolume> vol_fcn = new LSFindCellVolume("VolFcn", fine_patch_hierarchy);
            ls_fcn->setDataOnPatchHierarchy(ls_idx, ls_var, fine_patch_hierarchy, loop_time);

            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<ITC> ghost_cell_comps(1);
            ghost_cell_comps[0] = ITC(ls_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            HierarchyGhostCellInterpolation hier_ghost_cell;
            hier_ghost_cell.initializeOperatorState(ghost_cell_comps, fine_patch_hierarchy);
            hier_ghost_cell.fillData(loop_time);

            vol_fcn->updateVolumeAndArea(vol_idx, vol_var, IBTK::invalid_index, nullptr, ls_idx, ls_var, true);
            for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = fine_patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CellData<NDIM, int>> rank_data = patch->getPatchData(rank_idx);
                    rank_data->fillAll(SAMRAI_MPI::getRank());
                }
            }
        }

        {
            // Coarsened fine patch hierarchy
            Pointer<LSFindCellVolume> vol_fcn = new LSFindCellVolume("VolFcn", coarsened_fine_patch_hierarchy);
            ls_fcn->setDataOnPatchHierarchy(ls_idx, ls_var, coarsened_fine_patch_hierarchy, loop_time);

            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<ITC> ghost_cell_comps(1);
            ghost_cell_comps[0] = ITC(ls_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            HierarchyGhostCellInterpolation hier_ghost_cell;
            hier_ghost_cell.initializeOperatorState(ghost_cell_comps, coarsened_fine_patch_hierarchy);
            hier_ghost_cell.fillData(loop_time);

            vol_fcn->updateVolumeAndArea(vol_idx, vol_var, IBTK::invalid_index, nullptr, ls_idx, ls_var, true);
            for (int ln = 0; ln <= coarsened_fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CellData<NDIM, int>> rank_data = patch->getPatchData(rank_idx);
                    rank_data->fillAll(SAMRAI_MPI::getRank());
                }
            }
        }

        for (int ln = coarse_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            Pointer<PatchLevel<NDIM>> level = coarse_patch_hierarchy->getPatchLevel(ln - 1);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);
                pout << "Checking if data is valid on patch.\n";
            }
        }

        // Synchronize the coarse hierarchy data.
        for (int ln = coarse_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            Pointer<PatchLevel<NDIM>> coarser_level = coarse_patch_hierarchy->getPatchLevel(ln - 1);
            Pointer<PatchLevel<NDIM>> finer_level = coarse_patch_hierarchy->getPatchLevel(ln);

            CoarsenAlgorithm<NDIM> coarsen_alg;
            Pointer<CoarsenOperator<NDIM>> coarsen_op;

            coarsen_op = grid_geom->lookupCoarsenOperator(Q_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(Q_idx, Q_idx, coarsen_op);

            coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
        }

        // Synchronize the fine hierarchy data.
        for (int ln = fine_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            Pointer<PatchLevel<NDIM>> coarser_level = fine_patch_hierarchy->getPatchLevel(ln - 1);
            Pointer<PatchLevel<NDIM>> finer_level = fine_patch_hierarchy->getPatchLevel(ln);

            CoarsenAlgorithm<NDIM> coarsen_alg;
            Pointer<CoarsenOperator<NDIM>> coarsen_op;

            coarsen_op = grid_geom->lookupCoarsenOperator(Q_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(Q_idx, Q_idx, coarsen_op);

            coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
        }

        // Coarsen data from the fine hierarchy to the coarsened fine hierarchy.
        for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> dst_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM>> src_level = fine_patch_hierarchy->getPatchLevel(ln);

            Pointer<CoarsenOperator<NDIM>> coarsen_op;
            for (PatchLevel<NDIM>::Iterator p(dst_level); p; p++)
            {
                Pointer<Patch<NDIM>> dst_patch = dst_level->getPatch(p());
                Pointer<Patch<NDIM>> src_patch = src_level->getPatch(p());
                const Box<NDIM>& coarse_box = dst_patch->getBox();
                TBOX_ASSERT(Box<NDIM>::coarsen(src_patch->getBox(), 2) == coarse_box);

                coarsen_op = grid_geom->lookupCoarsenOperator(Q_var, "CONSERVATIVE_COARSEN");
                coarsen_op->coarsen(*dst_patch, *src_patch, Q_interp_idx, Q_idx, coarse_box, 2);
            }
        }

        // Interpolate and copy data from the coarsened fine patch hierarchy to
        // the coarse patch hierarchy.
        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            pout << "Interpolating on level " << ln << "\n";
            Pointer<PatchLevel<NDIM>> dst_level = coarse_patch_hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM>> src_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);

            RefineAlgorithm<NDIM> refine_alg;
            Pointer<RefineOperator<NDIM>> refine_op;

            refine_op = grid_geom->lookupRefineOperator(Q_var, "CONSERVATIVE_LINEAR_REFINE");
            refine_alg.registerRefine(Q_interp_idx, Q_interp_idx, Q_scr_idx, refine_op);

            ComponentSelector data_indices;
            data_indices.setFlag(Q_scr_idx);
            CartExtrapPhysBdryOp bc_helper(data_indices, "LINEAR");

            refine_alg.createSchedule(dst_level, src_level, ln - 1, coarsened_fine_patch_hierarchy, &bc_helper)
                ->fillData(loop_time);
        }

        // Output plot data before taking norms of differences.
        visit_data_writer->writePlotData(coarse_patch_hierarchy, 0, loop_time);
        visit_data_writer->writePlotData(coarsened_fine_patch_hierarchy, 1, loop_time);
        visit_data_writer->writePlotData(fine_patch_hierarchy, 2, loop_time);

        // Compute norms of differences.
        coarse_hier_cc_data_ops.subtract(Q_idx, Q_idx, Q_interp_idx);

        Pointer<LSFindCellVolume> vol_fcn = new LSFindCellVolume("VolFcn", coarse_patch_hierarchy);
        vol_fcn->updateVolumeAndArea(vol_idx, vol_var, IBTK::invalid_index, nullptr, ls_idx, ls_var, true);

        //        coarse_hier_cc_data_ops.multiply(wgt_cc_idx, vol_in_idx, wgt_cc_idx);
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
                    (*wgt_data)(idx) *= next_to_cut ? 0.0 : 1.0;
                }
            }
        }

        pout << "\n"
             << "Error in " << Q_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(Q_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(Q_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << coarse_hier_cc_data_ops.maxNorm(Q_idx, wgt_cc_idx) << "\n";

        // Output plot data after taking norms of differences.
        visit_data_writer->writePlotData(coarse_patch_hierarchy, 3, loop_time);

        pout << endl;
        pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout << endl;
    }

    SAMRAIManager::shutdown();
    SAMRAI_MPI::finalize();
    return 0;
} // main
