#ifndef included_CCAD_CutCellVolumeMeshMapping
#define included_CCAD_CutCellVolumeMeshMapping
#include "CCAD/CutCellMeshMapping.h"
#include "CCAD/FEMeshPartitioner.h"
#include "CCAD/ls_functions.h"
#include "CCAD/ls_utilities.h"

#include "ibtk/FEDataManager.h"

#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh.h"

namespace CCAD
{
/*!
 * CutCellVolumeMeshMapping maintains a description of the Lagrangian mesh from the point of view of the background
 * mesh. We maintain a mapping from each cut cell index to a vector of element and element parent pairs.
 */
class CutCellVolumeMeshMapping : public CutCellMeshMapping
{
public:
    /*!
     * \brief Constructor.
     */
    CutCellVolumeMeshMapping(std::string object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             const std::shared_ptr<FEMeshPartitioner>& fe_mesh_partitioner);

    /*!
     * \brief Constructor.
     */
    CutCellVolumeMeshMapping(std::string object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             const std::vector<std::shared_ptr<FEMeshPartitioner>>& fe_mesh_paritioners);

    /*!
     * \brief Constructor.
     */
    CutCellVolumeMeshMapping(std::string object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             IBTK::FEDataManager* fe_data_manager);

    /*!
     * \brief Constructor.
     */
    CutCellVolumeMeshMapping(std::string object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             const std::vector<IBTK::FEDataManager*>& fe_data_managers);

    /*!
     * \brief Default deconstructor.
     */
    virtual ~CutCellVolumeMeshMapping();

    /*!
     * \brief Deleted default constructor.
     */
    CutCellVolumeMeshMapping() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    CutCellVolumeMeshMapping(const CutCellVolumeMeshMapping& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    CutCellVolumeMeshMapping& operator=(const CutCellVolumeMeshMapping& that) = delete;

    void generateCutCellMappings() override;

    inline void registerMappingFcn(MappingFcn fcn, unsigned int part = 0)
    {
        d_mapping_fcns[part] = fcn;
    }

    inline std::shared_ptr<FEMeshPartitioner>& getMeshPartitioner(unsigned int part = 0)
    {
        return d_bdry_mesh_partitioners[part];
    }

private:
    void commonConstructor(const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db);

    bool findIntersection(libMesh::Point& p, libMesh::Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q);

    std::vector<std::shared_ptr<FEMeshPartitioner>> d_bdry_mesh_partitioners;
    std::vector<IBTK::FEDataManager*> d_fe_data_managers;
    unsigned int d_num_parts;

    std::vector<MappingFcn> d_mapping_fcns;
};

} // namespace CCAD
#endif
