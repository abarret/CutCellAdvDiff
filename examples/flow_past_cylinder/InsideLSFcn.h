// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2018 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_InsideLSFcn
#define included_InsideLSFcn

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/FESurfaceDistanceEvaluator.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>

#include <libmesh/boundary_mesh.h>

#include <CartesianGridGeometry.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the advected scalar Q.
 */
class InsideLSFcn : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    InsideLSFcn(const std::string& object_name,
                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                libMesh::BoundaryMesh* ls_mesh,
                libMesh::BoundaryMesh* disk_mesh,
                IBTK::FEDataManager* disk_fe_manager);

    /*!
     * \brief Destructor.
     */
    ~InsideLSFcn();

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool isTimeDependent() const override
    {
        return true;
    }

    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    /*!
     * Set the data on the patch interior to the exact answer.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(nullptr)) override;

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    InsideLSFcn();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    InsideLSFcn(const InsideLSFcn& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    InsideLSFcn& operator=(const InsideLSFcn& that);

    void updateMesh(libMesh::Mesh& ls_mesh, libMesh::Mesh& disk_mesh, IBTK::FEDataManager* vol_fe_manager);

    void initialLSValue(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM>> data,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time);

    int d_ls_scr_idx = IBTK::invalid_index, d_n_scr_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_ls_scr_var;

    libMesh::BoundaryMesh* d_ls_mesh;

    libMesh::BoundaryMesh* d_disk_mesh;
    IBTK::FEDataManager* d_disk_fe_manager;

    double d_cached_time = std::numeric_limits<double>::quiet_NaN();

    SAMRAI::tbox::Pointer<IBAMR::FESurfaceDistanceEvaluator> d_surface_distance_evaluator;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    IBTK::VectorNd d_cent;
    double d_R = std::numeric_limits<double>::quiet_NaN();

    static double s_large_num;
};
#endif //#ifndef included_InsideLSFcn
