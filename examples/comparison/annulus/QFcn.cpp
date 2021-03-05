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

#include "QFcn.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBAMR_config.h>

#include "LS/utility_functions.h"

#include <SAMRAI_config.h>

#include <array>

namespace LS
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

QFcn::QFcn(const string& object_name, Pointer<GridGeometry<NDIM>> grid_geom, Pointer<Database> input_db)
    : LSCartGridFunction(object_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif

    // Initialize object with data read from the input database.
    getFromInput(input_db);
    return;
} // QFcn

void
QFcn::setDataOnPatchHierarchy(const int data_idx,
                              Pointer<Variable<NDIM>> var,
                              Pointer<PatchHierarchy<NDIM>> hierarchy,
                              const double data_time,
                              const bool /*initial_time*/,
                              int coarsest_ln,
                              int finest_ln)
{
    coarsest_ln = coarsest_ln < 0 ? 0 : coarsest_ln;
    finest_ln = finest_ln < 0 ? hierarchy->getFinestLevelNumber() : finest_ln;

    auto integrator = IntegrateFunction::getIntegrator();

    auto fcn = [this](VectorNd X, double /*t*/) -> double {
        auto w = [](double r, double m) -> double {
            return std::cos(4.0 * r) + 3.0 / 4.0 * m * M_PI * std::log(r) -
                   3.0 / 32.0 * m * M_PI * (9.0 * std::log(3.0 * M_PI / 4.0) - std::log(M_PI / 4.0) - 4.0);
        };
        X = X - d_center;
        double r = X.norm();
        return w(r, d_m);
    };
    integrator->integrateFcnOnPatchHierarchy(hierarchy, d_ls_idx, data_idx, fcn, data_time);

    // Divide by total volume to get cell average
    for (int ln = coarsest_ln; ln <= finest_ln; ln++)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(data_idx);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
            Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                if ((*vol_data)(idx) > 0.0)
                {
                    VectorNd X = LS::find_cell_centroid(idx, *ls_data);
                    (*Q_data)(idx) /= (*vol_data)(idx)*dx[0] * dx[1]
#if (NDIM == 3)
                                      * dx[2]
#endif
                        ;
                }
                else
                {
                    (*Q_data)(idx) = 0.0;
                }
            }
        }
    }
}

void
QFcn::setDataOnPatch(const int data_idx,
                     Pointer<Variable<NDIM>> /*var*/,
                     Pointer<Patch<NDIM>> patch,
                     const double /*data_time*/,
                     const bool /*initial_time*/,
                     Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
#endif

    Q_data->fillAll(0.0);
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
QFcn::getFromInput(Pointer<Database> db)
{
    db->getDoubleArray("Center", d_center.data(), NDIM);
    d_m = db->getDouble("m");
    return;
} // getFromInput

} // namespace LS
