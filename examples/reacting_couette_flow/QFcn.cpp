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

#include "ibamr/app_namespaces.h"

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

    auto fcn = [this](VectorNd X, double t) -> double {
        auto w = [](double r, double D, double t) -> double {
            if (r < 1.0)
                return std::pow(std::cos(M_PI * r) + 1.0, 2.0);
            else
                return 0.0;
        };
        MatrixNd Q;
        Q(0, 0) = Q(1, 1) = std::cos(t * (2.0 * M_PI));
        Q(0, 1) = -std::sin(t * (2.0 * M_PI));
        Q(1, 0) = std::sin(t * (2.0 * M_PI));
        VectorNd current_com = Q * d_com;
        X = X - current_com;
        double r = X.norm();
        return w(r, d_D, t);
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
    d_D = db->getDouble("D");
    db->getDoubleArray("com", d_com.data(), NDIM);
    d_R = db->getDouble("r");
    return;
} // getFromInput

} // namespace LS
