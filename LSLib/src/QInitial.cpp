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

#include "LS/QInitial.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <SAMRAI_config.h>

#include <array>

namespace LS
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

QInitial::QInitial(const string& object_name, Pointer<GridGeometry<NDIM>> grid_geom, Pointer<Database> input_db)
    : LSCartGridFunction(object_name), d_grid_geom(grid_geom)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
    TBOX_ASSERT(d_grid_geom);
#endif

    // Initialize object with data read from the input database.
    getFromInput(input_db);
    return;
} // QInitial

QInitial::~QInitial()
{
    // intentionally blank
    return;
} // ~QInitial

void
QInitial::setDataOnPatchHierarchy(const int data_idx,
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

    if (d_init_type == "RADIAL")
    {
        auto fcn = [this](VectorNd X, double t) -> double {
            auto w = [](double r, double D, double t) -> double {
                return 10 * std::exp(-r * r / (4.0 * D * (t + 0.5))) / (4.0 * D * (t + 0.5));
            };
            X = X - d_center;
            for (int d = 0; d < NDIM; ++d) X(d) -= t * d_vel[d];
            double r = X.norm();
            return w(r, d_D, t);
        };
        integrator->integrateFcnOnPatchHierarchy(hierarchy, d_ls_idx, data_idx, fcn, data_time);
    }
    else if (d_init_type == "CONSTANT")
    {
        auto fcn = [this](VectorNd X, double t) -> double { return d_val; };
        integrator->integrateFcnOnPatchHierarchy(hierarchy, d_ls_idx, data_idx, fcn, data_time);
    }
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
                    (*Q_data)(idx) /= (*vol_data)(idx)*dx[0] * dx[1]
#if (NDIM == 3)
                                      * dx[2]
#endif
                        ;
                }
            }
        }
    }
}

void
QInitial::setDataOnPatch(const int data_idx,
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

void
QInitial::updateAverageToTotal(Pointer<CellVariable<NDIM, double>> Q_var,
                               Pointer<VariableContext> Q_ctx,
                               Pointer<PatchHierarchy<NDIM>> hierarchy,
                               const int vol_idx)
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, Q_ctx);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();

            Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(vol_idx);

            for (CellIterator<NDIM> ci(box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                (*Q_data)(idx) *= (*vol_data)(idx);
            }
        }
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
QInitial::getFromInput(Pointer<Database> db)
{
    if (db)
    {
        d_init_type = db->getStringWithDefault("init_type", d_init_type);

        if (d_init_type == "RADIAL")
        {
            d_R1 = db->getDoubleWithDefault("R1", d_R1);
            if (db->keyExists("Center")) db->getDoubleArray("Center", d_center.data(), NDIM);
            db->getDoubleArray("velocity", d_vel.data(), NDIM);
            d_D = db->getDoubleWithDefault("D", d_D);
        }
        else if (d_init_type == "ROTATIONAL")
        {
            d_R1 = db->getDoubleWithDefault("R1", d_R1);
            if (db->keyExists("Center")) db->getDoubleArray("Center", d_center.data(), NDIM);
            d_v = db->getDouble("v");
            d_D = db->getDoubleWithDefault("D", d_D);
        }
        else if (d_init_type == "CONSTANT")
        {
            d_val = db->getDouble("value");
        }
        else
        {
            TBOX_ERROR(d_object_name << "::getFromInput()\n"
                                     << "  invalid initialization type " << d_init_type << "\n");
        }
    }
    return;
} // getFromInput

} // namespace LS
