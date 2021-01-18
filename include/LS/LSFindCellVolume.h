#ifndef included_LS_LSFindCellVolume
#define included_LS_LSFindCellVolume

#include "CellVariable.h"
#include "NodeVariable.h"
#include "PatchHierarchy.h"
#include "SideVariable.h"
#include "tbox/Pointer.h"

namespace LS
{
class LSFindCellVolume : public virtual SAMRAI::tbox::DescribedClass
{
public:
    LSFindCellVolume(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    virtual ~LSFindCellVolume() = default;

    virtual void updateVolumeAreaSideLS(int vol_idx,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vol_var,
                                        int area_idx,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> area_var,
                                        int side_idx,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> side_var,
                                        int phi_idx,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> phi_var,
                                        double data_time,
                                        bool extended_box = false) = 0;

    inline void setLS(bool set_ls)
    {
        d_set_ls = set_ls;
    }

protected:
    std::string d_object_name;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    bool d_set_ls = true;
};
} // namespace LS

#endif
