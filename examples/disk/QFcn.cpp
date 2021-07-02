#include "ibamr/config.h"

#include "CCAD/app_namespaces.h"

#include "QFcn.h"

#include <SAMRAI_config.h>

#include <array>

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

    auto fcn = [this](VectorNd X, double t) -> double {
        X -= d_cent;
        return 1.0 + X.squaredNorm() * (-d_sf_max * d_k_on + (d_k_off + d_k_on) * t * (1.0 + t)) /
                         (2.0 * d_D + d_k_on * (d_sf_max - t * (1.0 + t)));
    };

    // Divide by total volume to get cell average
    for (int ln = coarsest_ln; ln <= finest_ln; ln++)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(data_idx);
            Q_data->fillAll(0.0);
            Pointer<CellData<NDIM, double>> vol_data = patch->getPatchData(d_vol_idx);
            Pointer<NodeData<NDIM, double>> ls_data = patch->getPatchData(d_ls_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                if ((*vol_data)(idx) > 0.0)
                {
                    VectorNd cell_centroid = CCAD::find_cell_centroid(idx, *ls_data);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        cell_centroid[d] = xlow[d] + dx[d] * cell_centroid[d];
                    }
                    (*Q_data)(idx) = fcn(cell_centroid, data_time);
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
    d_D = db->getDouble("d");
    d_k_off = db->getDouble("k_off");
    d_k_on = db->getDouble("k_on");
    d_sf_max = db->getDouble("sf_max");
    db->getDoubleArray("center", d_cent.data(), NDIM);
    return;
} // getFromInput
