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

#ifndef included_ForcingFcn
#define included_ForcingFcn

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

#include <CartesianGridGeometry.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include "LS/IntegrateFunction.h"
#include "LS/LSCartGridFunction.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the advected scalar Q.
 */
namespace LS
{
class ForcingFcn : public LS::LSCartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    ForcingFcn(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~ForcingFcn();

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool isTimeDependent() const
    {
        return true;
    }

    /*!
     * Set the data on the patch interior to the exact answer. Note that we are setting the cell average here.
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
    ForcingFcn();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    ForcingFcn(const ForcingFcn& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ForcingFcn& operator=(const ForcingFcn& that);

    double d_D = 0.0;
    double d_k_on = 0.0;
    double d_k_off = 0.0;
    double d_sf_max = 0.0;
    IBTK::VectorNd d_cent = { 0.0, 0.0 };
};

} // namespace LS

#endif //#ifndef included_ForcingFcn
