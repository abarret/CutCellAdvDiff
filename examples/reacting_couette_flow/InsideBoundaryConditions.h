#ifndef included_InsideBoundaryConditions
#define included_InsideBoundaryConditions

#include "LS/LSCutCellBoundaryConditions.h"
#include "LS/SemiLagrangianAdvIntegrator.h"

class InsideBoundaryConditions : public LS::LSCutCellBoundaryConditions
{
public:
    InsideBoundaryConditions(const std::string& object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> out_var,
                             SAMRAI::tbox::Pointer<LS::SemiLagrangianAdvIntegrator> integrator);

    ~InsideBoundaryConditions();

    inline void setContext(SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
    {
        d_ctx = ctx;
    }

    /*!
     * \brief Deleted default constructor.
     */
    InsideBoundaryConditions() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    InsideBoundaryConditions(const InsideBoundaryConditions& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    InsideBoundaryConditions& operator=(const InsideBoundaryConditions& that) = delete;

    void applyBoundaryCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                int Q_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> R_var,
                                int R_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                double time);

private:
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_out_var;
    SAMRAI::tbox::Pointer<LS::SemiLagrangianAdvIntegrator> d_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;

    double d_k1 = std::numeric_limits<double>::quiet_NaN();
    double d_D_coef = std::numeric_limits<double>::quiet_NaN();
};
#endif
