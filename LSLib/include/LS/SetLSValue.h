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

#ifndef included_LS_SetLSValue
#define included_LS_SetLSValue

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

// C++ namespace delcarations
#include <ibamr/app_namespaces.h>

namespace LS
{
/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the advected scalar Q.
 */
class SetLSValue : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    SetLSValue(const string& object_name, Pointer<GridGeometry<NDIM>> grid_geom, Pointer<Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~SetLSValue();

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
    SetLSValue();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SetLSValue(const SetLSValue& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SetLSValue& operator=(const SetLSValue& that);

    /*!
     * Read input values, indicated above, from given database.
     */
    void getFromInput(Pointer<Database> db);

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    string d_object_name;

    /*
     * The grid geometry.
     */
    Pointer<CartesianGridGeometry<NDIM>> d_grid_geom;

    /*
     * The initialization type.
     */
    string d_interface_type = "ANNULUS";

    IBTK::VectorNd d_U;

    /*
     * Disk information
     */
    double d_R1 = 0.25;
#if (NDIM == 2)
    VectorNd d_center = { 1.509, 1.521 };
#endif
#if (NDIM == 3)
    VectorNd d_center = { 1.509, 1.521, 1.514 };
#endif

    bool d_extended_box = false;
};
} // namespace LS

#endif //#ifndef included_LS_SetLSValue
