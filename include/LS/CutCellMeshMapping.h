#ifndef included_CutCellMeshMapping
#define included_CutCellMeshMapping
#include "ibtk/FEDataManager.h"

#include "LS/utility_functions.h"

#include "libmesh/mesh.h"

namespace LS
{
/*!
 * CutCellMeshMapping maintains a description of the Lagrangian mesh from a description of the background mesh. We
 * maintain a mapping from each cut cell index to a vector of element and element parent pairs.
 */
class CutCellMeshMapping
{
public:
    /*!
     * \brief Constructor.
     */
    CutCellMeshMapping(std::string object_name,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                       libMesh::Mesh* mesh,
                       IBTK::FEDataManager* fe_data_manager);

    /*!
     * \brief Default deconstructor.
     */
    ~CutCellMeshMapping() = default;

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

    void setLSData(int ls_idx, int vol_idx, int area_idx);

    void initializeObjectState(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    void deinitializeObjectState();

    void generateCutCellMappings();

    inline const std::map<LS::PatchIndexPair, std::vector<CutCellElems>>& getIdxCutCellElemsMap(const int ln)
    {
        return d_idx_cut_cell_elems_map_vec[ln];
    }

private:
    bool findIntersection(libMesh::Point& p, libMesh::Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q);

    std::string d_object_name;
    bool d_is_initialized = false;
    bool d_perturb_nodes = false;

    libMesh::Mesh* d_mesh = nullptr;
    IBTK::FEDataManager* d_fe_data_manager = nullptr;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    int d_ls_idx = IBTK::invalid_index, d_vol_idx = IBTK::invalid_index, d_area_idx = IBTK::invalid_index;

    std::vector<std::map<LS::PatchIndexPair, std::vector<CutCellElems>>> d_idx_cut_cell_elems_map_vec;
};

} // namespace LS
#endif
