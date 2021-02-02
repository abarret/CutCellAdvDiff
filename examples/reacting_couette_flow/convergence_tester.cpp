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
#include <IBAMR_config.h>
#include <IBTK_config.h>

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
#include "LS/utility_functions.h"

#include "BoxArray.h"
#include "CartesianPatchGeometry.h"
#include "CoarseFineBoundary.h"
#include "InsideLSFcn.h"
#include "OutsideLSFcn.h"
#include "PatchGeometry.h"
#include "RefineAlgorithm.h"

using namespace LS;

static double a = std::numeric_limits<double>::signaling_NaN();
static double b = std::numeric_limits<double>::signaling_NaN();

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

    a = input_db->getDouble("A");
    b = input_db->getDouble("B");

    // Create major algorithm and data objects that comprise application.
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    // Initialize variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    Pointer<VariableContext> current_ctx = var_db->getContext("SemiLagrangianAdvIntegrator::CURRENT");
    Pointer<VariableContext> scratch_ctx = var_db->getContext("SemiLagrangianAdvIntegrator::SCRATCH");

    Pointer<CellVariable<NDIM, double>> Q_in_var = new CellVariable<NDIM, double>("Q_in");
    const int Q_in_idx = var_db->registerVariableAndContext(Q_in_var, current_ctx);
    const int Q_in_interp_idx = var_db->registerClonedPatchDataIndex(Q_in_var, Q_in_idx);
    const int Q_in_scr_idx = var_db->registerVariableAndContext(Q_in_var, scratch_ctx, 2);

    Pointer<CellVariable<NDIM, double>> Q_out_var = new CellVariable<NDIM, double>("Q_out");
    const int Q_out_idx = var_db->registerVariableAndContext(Q_out_var, current_ctx);
    const int Q_out_interp_idx = var_db->registerClonedPatchDataIndex(Q_out_var, Q_out_idx);
    const int Q_out_scr_idx = var_db->registerVariableAndContext(Q_out_var, scratch_ctx, 2);

    Pointer<NodeVariable<NDIM, double>> ls_in_var = new NodeVariable<NDIM, double>("LS_in_var");
    const int ls_in_idx = var_db->registerVariableAndContext(ls_in_var, current_ctx, IntVector<NDIM>(2));

    Pointer<NodeVariable<NDIM, double>> ls_out_var = new NodeVariable<NDIM, double>("LS_out_var");
    const int ls_out_idx = var_db->registerVariableAndContext(ls_out_var, current_ctx, IntVector<NDIM>(2));

    Pointer<CellVariable<NDIM, double>> vol_in_var = new CellVariable<NDIM, double>("Volume In"),
                                        vol_out_var = new CellVariable<NDIM, double>("Volume Out");
    const int vol_in_idx = var_db->registerVariableAndContext(vol_in_var, current_ctx, IntVector<NDIM>(2));
    const int vol_out_idx = var_db->registerVariableAndContext(vol_out_var, current_ctx, IntVector<NDIM>(2));

    Pointer<CellVariable<NDIM, double>> area_in_var = new CellVariable<NDIM, double>("Area In"),
                                        area_out_var = new CellVariable<NDIM, double>("Area Out");
    const int area_in_idx = var_db->registerVariableAndContext(area_in_var, current_ctx, IntVector<NDIM>(2));
    const int area_out_idx = var_db->registerVariableAndContext(area_out_var, current_ctx, IntVector<NDIM>(2));

    Pointer<CellVariable<NDIM, int>> rank_var = new CellVariable<NDIM, int>("Rank");
    const int rank_idx = var_db->registerVariableAndContext(rank_var, current_ctx);

    // Set up visualization plot file writer.
    Pointer<VisItDataWriter<NDIM>> visit_data_writer =
        new VisItDataWriter<NDIM>("VisIt Writer", main_db->getString("viz_dump_dirname"), 1);
    visit_data_writer->registerPlotQuantity("Q_in", "SCALAR", Q_in_idx);
    visit_data_writer->registerPlotQuantity("Q_in interp", "SCALAR", Q_in_interp_idx);
    visit_data_writer->registerPlotQuantity("Q_out", "SCALAR", Q_out_idx);
    visit_data_writer->registerPlotQuantity("Q_out interp", "SCALAR", Q_out_interp_idx);
    visit_data_writer->registerPlotQuantity("ls_in", "SCALAR", ls_in_idx);
    visit_data_writer->registerPlotQuantity("ls_out", "SCALAR", ls_out_idx);
    visit_data_writer->registerPlotQuantity("vol_in", "SCALAR", vol_in_idx);
    visit_data_writer->registerPlotQuantity("vol_out", "SCALAR", vol_out_idx);
    visit_data_writer->registerPlotQuantity("area_in", "SCALAR", area_in_idx);
    visit_data_writer->registerPlotQuantity("area_out", "SCALAR", area_out_idx);
    visit_data_writer->registerPlotQuantity("rank", "SCALAR", rank_idx);

    // Time step loop.
    double loop_time = 0.0;
    int coarse_iteration_num = coarse_hier_dump_interval;
    int fine_iteration_num = fine_hier_dump_interval;

    bool files_exist = true;
    int draw_iter = 0;
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
        hier_data.setFlag(Q_in_idx);
        hier_data.setFlag(Q_out_idx);

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
            level->allocatePatchData(Q_in_interp_idx, loop_time);
            level->allocatePatchData(Q_out_interp_idx, loop_time);
            level->allocatePatchData(ls_in_idx, loop_time);
            level->allocatePatchData(ls_out_idx, loop_time);
            level->allocatePatchData(vol_in_idx, loop_time);
            level->allocatePatchData(vol_out_idx, loop_time);
            level->allocatePatchData(area_in_idx, loop_time);
            level->allocatePatchData(area_out_idx, loop_time);
            level->allocatePatchData(Q_in_scr_idx, loop_time);
            level->allocatePatchData(Q_out_scr_idx, loop_time);
            level->allocatePatchData(rank_idx, loop_time);
        }

        for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = fine_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_in_interp_idx, loop_time);
            level->allocatePatchData(Q_out_interp_idx, loop_time);
            level->allocatePatchData(Q_in_scr_idx, loop_time);
            level->allocatePatchData(Q_out_scr_idx, loop_time);
            level->allocatePatchData(ls_in_idx, loop_time);
            level->allocatePatchData(ls_out_idx, loop_time);
            level->allocatePatchData(vol_in_idx, loop_time);
            level->allocatePatchData(vol_out_idx, loop_time);
            level->allocatePatchData(area_in_idx, loop_time);
            level->allocatePatchData(area_out_idx, loop_time);
            level->allocatePatchData(rank_idx, loop_time);
        }

        for (int ln = 0; ln <= coarsened_fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_in_interp_idx, loop_time);
            level->allocatePatchData(Q_out_interp_idx, loop_time);
            level->allocatePatchData(Q_in_idx, loop_time);
            level->allocatePatchData(Q_out_idx, loop_time);
            level->allocatePatchData(Q_in_scr_idx, loop_time);
            level->allocatePatchData(Q_out_scr_idx, loop_time);
            level->allocatePatchData(ls_in_idx, loop_time);
            level->allocatePatchData(ls_out_idx, loop_time);
            level->allocatePatchData(vol_in_idx, loop_time);
            level->allocatePatchData(vol_out_idx, loop_time);
            level->allocatePatchData(rank_idx, loop_time);
        }

        Pointer<CartGridFunction> ls_in_fcn =
            new InsideLSFcn("InsideLSFcn", app_initializer->getComponentDatabase("InsideLSFcn"));
        Pointer<CartGridFunction> ls_out_fcn = new OutsideLSFcn(
            current_ctx, "OutsideLSFcn", nullptr, ls_in_var, app_initializer->getComponentDatabase("OutsideLSFcn"));
        {
            // Coarse patch hierarchy
            Pointer<LSFindCellVolume> vol_fcn = new LSFindCellVolume("VolFcn", coarse_patch_hierarchy);
            ls_in_fcn->setDataOnPatchHierarchy(ls_in_idx, ls_in_var, coarse_patch_hierarchy, loop_time);
            ls_out_fcn->setDataOnPatchHierarchy(ls_out_idx, ls_out_var, coarse_patch_hierarchy, loop_time);

            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<ITC> ghost_cell_comps(2);
            ghost_cell_comps[0] = ITC(ls_in_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            ghost_cell_comps[1] = ITC(ls_out_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            HierarchyGhostCellInterpolation hier_ghost_cell;
            hier_ghost_cell.initializeOperatorState(ghost_cell_comps, coarse_patch_hierarchy);
            hier_ghost_cell.fillData(loop_time);

            vol_fcn->updateVolumeAndArea(vol_in_idx, vol_in_var, area_in_idx, area_in_var, ls_in_idx, ls_in_var, true);
            vol_fcn->updateVolumeAndArea(
                vol_out_idx, vol_out_var, area_out_idx, area_out_var, ls_out_idx, ls_out_var, true);
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
            ls_in_fcn->setDataOnPatchHierarchy(ls_in_idx, ls_in_var, fine_patch_hierarchy, loop_time);
            ls_out_fcn->setDataOnPatchHierarchy(ls_out_idx, ls_out_var, fine_patch_hierarchy, loop_time);

            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<ITC> ghost_cell_comps(2);
            ghost_cell_comps[0] = ITC(ls_in_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            ghost_cell_comps[1] = ITC(ls_out_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            HierarchyGhostCellInterpolation hier_ghost_cell;
            hier_ghost_cell.initializeOperatorState(ghost_cell_comps, fine_patch_hierarchy);
            hier_ghost_cell.fillData(loop_time);

            vol_fcn->updateVolumeAndArea(vol_in_idx, vol_in_var, area_in_idx, area_in_var, ls_in_idx, ls_in_var, true);
            vol_fcn->updateVolumeAndArea(
                vol_out_idx, vol_out_var, area_out_idx, area_out_var, ls_out_idx, ls_out_var, true);
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
            ls_in_fcn->setDataOnPatchHierarchy(ls_in_idx, ls_in_var, coarsened_fine_patch_hierarchy, loop_time);
            ls_out_fcn->setDataOnPatchHierarchy(ls_out_idx, ls_out_var, coarsened_fine_patch_hierarchy, loop_time);

            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<ITC> ghost_cell_comps(2);
            ghost_cell_comps[0] = ITC(ls_in_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            ghost_cell_comps[1] = ITC(ls_out_idx, "LINEAR_REFINE", false, "NONE", "LINEAR");
            HierarchyGhostCellInterpolation hier_ghost_cell;
            hier_ghost_cell.initializeOperatorState(ghost_cell_comps, coarsened_fine_patch_hierarchy);
            hier_ghost_cell.fillData(loop_time);

            vol_fcn->updateVolumeAndArea(
                vol_in_idx, vol_in_var, IBTK::invalid_index, nullptr, ls_in_idx, ls_in_var, true);
            vol_fcn->updateVolumeAndArea(
                vol_out_idx, vol_out_var, IBTK::invalid_index, nullptr, ls_out_idx, ls_out_var, true);
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
                Pointer<CellData<NDIM, double>> vol_in_data = patch->getPatchData(vol_in_idx);
                Pointer<CellData<NDIM, double>> vol_out_data = patch->getPatchData(vol_out_idx);
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

            coarsen_op = grid_geom->lookupCoarsenOperator(Q_in_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(Q_in_idx, Q_in_idx, coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(Q_out_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(Q_out_idx, Q_out_idx, coarsen_op);

            coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
        }

        // Synchronize the fine hierarchy data.
        for (int ln = fine_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            Pointer<PatchLevel<NDIM>> coarser_level = fine_patch_hierarchy->getPatchLevel(ln - 1);
            Pointer<PatchLevel<NDIM>> finer_level = fine_patch_hierarchy->getPatchLevel(ln);

            CoarsenAlgorithm<NDIM> coarsen_alg;
            Pointer<CoarsenOperator<NDIM>> coarsen_op;

            coarsen_op = grid_geom->lookupCoarsenOperator(Q_in_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(Q_in_idx, Q_in_idx, coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(Q_out_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(Q_out_idx, Q_out_idx, coarsen_op);

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

                coarsen_op = grid_geom->lookupCoarsenOperator(Q_in_var, "CONSERVATIVE_COARSEN");
                coarsen_op->coarsen(*dst_patch, *src_patch, Q_in_interp_idx, Q_in_idx, coarse_box, 2);

                coarsen_op = grid_geom->lookupCoarsenOperator(Q_out_var, "CONSERVATIVE_COARSEN");
                coarsen_op->coarsen(*dst_patch, *src_patch, Q_out_interp_idx, Q_out_idx, coarse_box, 2);
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

            refine_op = grid_geom->lookupRefineOperator(Q_in_var, "CONSERVATIVE_LINEAR_REFINE");
            refine_alg.registerRefine(Q_in_interp_idx, Q_in_interp_idx, Q_in_scr_idx, refine_op);

            refine_op = grid_geom->lookupRefineOperator(Q_out_var, "CONSERVATIVE_LINEAR_REFINE");
            refine_alg.registerRefine(Q_out_interp_idx, Q_out_interp_idx, Q_out_scr_idx, refine_op);

            ComponentSelector data_indices;
            data_indices.setFlag(Q_in_scr_idx);
            data_indices.setFlag(Q_out_scr_idx);
            CartExtrapPhysBdryOp bc_helper(data_indices, "LINEAR");

            refine_alg.createSchedule(dst_level, src_level, ln - 1, coarsened_fine_patch_hierarchy, &bc_helper)
                ->fillData(loop_time);
        }

        // Output boundary information
        outputBdryInfo(Q_in_idx,
                       Q_in_scr_idx,
                       ls_in_idx,
                       vol_in_idx,
                       area_in_idx,
                       loop_time,
                       coarse_iteration_num,
                       coarse_patch_hierarchy,
                       "coarse/in_bdry_info_");
        outputBdryInfo(Q_in_idx,
                       Q_in_scr_idx,
                       ls_in_idx,
                       vol_in_idx,
                       area_in_idx,
                       loop_time,
                       fine_iteration_num,
                       fine_patch_hierarchy,
                       "fine/in_bdry_info_");
        outputBdryInfo(Q_out_idx,
                       Q_out_scr_idx,
                       ls_out_idx,
                       vol_out_idx,
                       area_in_idx,
                       loop_time,
                       coarse_iteration_num,
                       coarse_patch_hierarchy,
                       "coarse/out_bdry_info_");
        outputBdryInfo(Q_out_idx,
                       Q_out_scr_idx,
                       ls_out_idx,
                       vol_out_idx,
                       area_in_idx,
                       loop_time,
                       fine_iteration_num,
                       fine_patch_hierarchy,
                       "fine/out_bdry_info_");

        // Output plot data before taking norms of differences.
        visit_data_writer->writePlotData(coarse_patch_hierarchy, draw_iter++, loop_time);
        visit_data_writer->writePlotData(coarsened_fine_patch_hierarchy, draw_iter++, loop_time);
        visit_data_writer->writePlotData(fine_patch_hierarchy, draw_iter++, loop_time);

        // Compute norms of differences.
        coarse_hier_cc_data_ops.subtract(Q_in_idx, Q_in_idx, Q_in_interp_idx);
        coarse_hier_cc_data_ops.subtract(Q_out_idx, Q_out_idx, Q_out_interp_idx);

        //        coarse_hier_cc_data_ops.multiply(wgt_cc_idx, vol_in_idx, wgt_cc_idx);
        const IntVector<NDIM> cut_cell_width(input_db->getInteger("SKIP_CUT_CELLS"));
        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = coarse_patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);
                Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_in_idx);
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

        pout << "\n"
             << "Error in " << Q_in_var->getName() << " at time " << loop_time << ":\n"
             << "  In_L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(Q_in_idx, wgt_cc_idx) << "\n"
             << "  In_L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(Q_in_idx, wgt_cc_idx) << "\n"
             << "  In_max-norm: " << coarse_hier_cc_data_ops.maxNorm(Q_in_idx, wgt_cc_idx) << "\n";

        hier_math_ops.resetLevels(0, coarse_patch_hierarchy->getFinestLevelNumber());
        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = coarse_patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);
                Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_out_idx);
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
        //        coarse_hier_cc_data_ops.multiply(wgt_cc_idx, vol_out_idx, wgt_cc_idx);

        pout << "\n"
             << "Error in " << Q_out_var->getName() << " at time " << loop_time << ":\n"
             << "  Out_L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(Q_out_idx, wgt_cc_idx) << "\n"
             << "  Out_L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(Q_out_idx, wgt_cc_idx) << "\n"
             << "  Out_max-norm: " << coarse_hier_cc_data_ops.maxNorm(Q_out_idx, wgt_cc_idx) << "\n";

        // Output plot data after taking norms of differences.
        visit_data_writer->writePlotData(coarse_patch_hierarchy, draw_iter++, loop_time);

        pout << endl;
        pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout << endl;
    }

    SAMRAIManager::shutdown();
    SAMRAI_MPI::finalize();
    return 0;
} // main

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
                    cur_theta + (a * cur_r + b / cur_r) * (2.0 / (2.0 * M_PI)) * std::cos(2 * M_PI * loop_time / 2.0) -
                    (a * cur_r + b / cur_r) * (2.0 / (2.0 * M_PI));
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
