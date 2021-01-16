#ifndef included_OutsideBoundaryConditions
#define included_OutsideBoundaryConditions

#include "LS/LSCutCellBoundaryConditions.h"
#include "LS/SemiLagrangianAdvIntegrator.h"

class OutsideBoundaryConditions : public LS::LSCutCellBoundaryConditions
{
public:
    OutsideBoundaryConditions(const std::string& object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> in_var,
                              SAMRAI::tbox::Pointer<LS::SemiLagrangianAdvIntegrator> integrator);

    ~OutsideBoundaryConditions();

    inline void setContext(SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
    {
        d_ctx = ctx;
    }

    /*!
     * \brief Deleted default constructor.
     */
    OutsideBoundaryConditions() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    OutsideBoundaryConditions(const OutsideBoundaryConditions& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    OutsideBoundaryConditions& operator=(const OutsideBoundaryConditions& that) = delete;

    void applyBoundaryCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                int Q_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> R_var,
                                int R_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                double time);

    inline void registerAreaAndLSInsideIndex(const int area_in_idx, int ls_in_idx)
    {
        d_ls_in_idx = ls_in_idx;
        d_area_in_idx = area_in_idx;
    }

private:
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_in_var;
    SAMRAI::tbox::Pointer<LS::SemiLagrangianAdvIntegrator> d_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;

    int d_ls_in_idx = IBTK::invalid_index;
    int d_area_in_idx = IBTK::invalid_index;

    double d_k1 = std::numeric_limits<double>::quiet_NaN();
};
#endif
