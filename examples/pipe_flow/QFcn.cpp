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

#include <ibamr/config.h>

#include "QFcn.h"

#include <SAMRAI_config.h>

#include <array>

namespace LS
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

QFcn::QFcn(const string& object_name, Pointer<Database> input_db) : LSCartGridFunction(object_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif

    // Initialize object with data read from the input database.
    getFromInput(input_db);
    return;
} // QFcn

QFcn::~QFcn()
{
    // intentionally blank
    return;
} // ~QFcn

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

    auto fcn = [this](VectorNd X, double t) -> double { return 0.0; };
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
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                if ((*vol_data)(idx) > 0.0)
                {
                    (*Q_data)(idx) /= (*vol_data)(idx)*dx[0] * dx[1];
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
    d_theta = db->getDouble("theta");
    db->getDoubleArray("channel_center", d_channel_center.data(), NDIM);
    d_y_low = db->getDouble("y_low");
    d_y_up = db->getDouble("y_up");
    d_D = db->getDouble("d");
    return;
} // getFromInput

} // namespace LS
