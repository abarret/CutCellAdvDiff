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

#ifndef included_OutsideLSFcn
#define included_OutsideLSFcn

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>

#include <CartesianGridGeometry.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the advected scalar Q.
 */
class OutsideLSFcn : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    OutsideLSFcn(const std::string& object_name,
                 SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> hierarchy_integrator,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> in_var,
                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    OutsideLSFcn(SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
                 const std::string& object_name,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> in_var,
                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~OutsideLSFcn();

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool isTimeDependent() const override
    {
        return true;
    }

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
    OutsideLSFcn();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    OutsideLSFcn(const OutsideLSFcn& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    OutsideLSFcn& operator=(const OutsideLSFcn& that);

    SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> d_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_in_ls_var;
    double d_R1 = std::numeric_limits<double>::quiet_NaN(), d_R2 = std::numeric_limits<double>::quiet_NaN();
};
#endif //#ifndef included_OutsideLSFcn
