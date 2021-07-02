#ifndef included_CCAD_CutCellMeshMapping
#define included_CCAD_CutCellMeshMapping
#include "CCAD/ls_functions.h"
#include "CCAD/ls_utilities.h"

#include "ibtk/FEDataManager.h"

#include "libmesh/mesh.h"

namespace CCAD
{
/*!
 * CutCellMeshMapping maintains a description of the Lagrangian mesh from the point of view of the background mesh. We
 * maintain a mapping from each cut cell index to a vector of element and element parent pairs.
 */
class CutCellMeshMapping : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    CutCellMeshMapping(std::string object_name,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                       unsigned int parts = 0);

    /*!
     * \brief Default deconstructor.
     */
    virtual ~CutCellMeshMapping();

    /*!
     * \brief Deleted default constructor.
     */
    CutCellMeshMapping() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    CutCellMeshMapping(const CutCellMeshMapping& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    CutCellMeshMapping& operator=(const CutCellMeshMapping& that) = delete;

    virtual void initializeObjectState(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    virtual void deinitializeObjectState();

    virtual void generateCutCellMappings() = 0;

    inline const std::vector<std::map<IndexList, std::vector<CutCellElems>>>& getIdxCutCellElemsMap(const int ln)
    {
        return d_idx_cut_cell_elems_map_vec[ln];
    }

    inline unsigned int getNumParts()
    {
        return d_num_parts;
    }

protected:
    std::string d_object_name;
    bool d_is_initialized = false;
    bool d_perturb_nodes = false;
    unsigned int d_num_parts = 0;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    std::vector<std::vector<std::map<IndexList, std::vector<CutCellElems>>>> d_idx_cut_cell_elems_map_vec;
};

} // namespace CCAD
#endif
