#ifndef included_SBSurfaceFluidCouplingManager
#define included_SBSurfaceFluidCouplingManager

#include "ibtk/FEDataManager.h"

#include "LS/CutCellMeshMapping.h"
#include "LS/RBFReconstructCache.h"

#include "CellData.h"
#include "NodeData.h"

#include "libmesh/equation_systems.h"
#include <libmesh/mesh.h>

namespace LS
{
class SBSurfaceFluidCouplingManager
{
public:
    /*!
     * \brief Constructor.
     */
    SBSurfaceFluidCouplingManager(std::string name,
                                  const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
                                  IBTK::FEDataManager* fe_data_manager,
                                  libMesh::Mesh* mesh);

    /*!
     * \brief Deconstructor. Cleans up any allocated data objects.
     */
    ~SBSurfaceFluidCouplingManager();

    /*!
     * \brief Deleted default constructor.
     */
    SBSurfaceFluidCouplingManager() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    SBSurfaceFluidCouplingManager(const SBSurfaceFluidCouplingManager& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    SBSurfaceFluidCouplingManager& operator=(const SBSurfaceFluidCouplingManager& that) = delete;

    /*!
     * \brief Register surface concentrations that will be tracked by the data manager.
     * @{
     */
    void registerSurfaceConcentration(std::string surface_name);
    void registerSurfaceConcentration(const std::vector<std::string>& surface_names);
    /*!
     * @}
     */

    /*!
     * \brief Register fluid phase concentrations. The data manager will keep a reference to the variable and provide an
     * interface to calculate concentrations along the surface.
     *
     * \see interpolateToBoundary
     * @{
     */
    void registerFluidConcentration(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> fl_var);
    void registerFluidConcentration(
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>>& fl_vars);
    /*!
     * @}
     */

    /*!
     * \brief Register dependence between the given surface concentration and fluid variable.
     */
    void registerFluidSurfaceDependence(const std::string& surface_name,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> fl_var);

    /*!
     * \brief Register dependence between two surface concentrations.
     */
    void registerSurfaceSurfaceDependence(const std::string& part1_name, const std::string& part2_name);

    /*!
     * \brief Register a reaction function and context for a given surface quantity.
     */
    void registerSurfaceReactionFunction(const std::string& surface_name, ReactionFcn fcn, void* ctx = nullptr);

    /*!
     * \brief Register boundary conditions for a given fluid quantity.
     */
    void registerFluidBoundaryCondition(const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>& fl_var,
                                        ReactionFcn a_fcn,
                                        ReactionFcn g_fcn,
                                        void* ctx = nullptr);

    /*!
     * \brief Initialize the equation system holding the surface and fluid concentrations. Note that all concentrations
     * should be registered before this function is called.
     *
     * Note: This function calls EquationSystems::reinit(), so any call to EquationSystems::init() MUST occur prior to
     * this function call.
     */
    void initializeFEEquationSystems();

    /*!
     * \brief Set the level set and volume data used by this manager to integrate fluid concentrations to the boundary.
     */
    void setLSData(int ls_idx, int vol_idx, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    /*!
     * \brief Interpolate the data stored in the SAMRAI variable and context pair to the surface.
     * @{
     */
    const std::string& interpolateToBoundary(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> fl_var,
                                             const int idx,
                                             double time);
    inline const std::string&
    interpolateToBoundary(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> fl_var,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
                          double time)
    {
        auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int idx = var_db->mapVariableAndContextToIndex(fl_var, ctx);
        return interpolateToBoundary(fl_var, idx, time);
    }
    inline const std::string& interpolateToBoundary(const std::string& fl_name,
                                                    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
                                                    double time)
    {
        auto fl_it = std::find(d_fl_names.begin(), d_fl_names.end(), fl_name);
        TBOX_ASSERT(fl_it != d_fl_names.end());
        return interpolateToBoundary(d_fl_vars[std::distance(d_fl_names.begin(), fl_it)], ctx, time);
    }
    inline const std::string& interpolateToBoundary(const std::string& fl_name, const int idx, double time)
    {
        auto fl_it = std::find(d_fl_names.begin(), d_fl_names.end(), fl_name);
        TBOX_ASSERT(fl_it != d_fl_names.end());
        return interpolateToBoundary(d_fl_vars[std::distance(d_fl_names.begin(), fl_it)], idx, time);
    }
    /*!
     * @}
     */

    /*!
     * \brief Update the Jacobian of the mapping. Returns the string of the system name.
     */
    const std::string& updateJacobian();

    /*!
     * \brief Return the string containing the Jacobian of the mapping.
     */
    inline const std::string& getJacobianName()
    {
        return d_J_sys_name;
    }

    /*!
     * \brief Returns the FEDataManager object used by this manager.
     */
    inline IBTK::FEDataManager* getFEDataManager()
    {
        return d_fe_data_manager;
    }

    /*!
     * \brief Return the surface concentration variable names stored by this manager.
     */
    inline const std::vector<std::string>& getSFNames()
    {
        return d_sf_names;
    }

    /*!
     * \brief Return the fluid concentration variable names stored by this manager.
     */
    inline const std::vector<std::string>& getFLNames()
    {
        return d_fl_names;
    }

    /*!
     * \brief Return the SAMRAI variable represented by the fluid variable name.
     */
    inline const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>&
    getFLVariable(const std::string& fl_name)
    {
        auto fl_it = std::find(d_fl_names.begin(), d_fl_names.end(), fl_name);
        TBOX_ASSERT(fl_it != d_fl_names.end());
        return d_fl_vars[std::distance(d_fl_names.begin(), fl_it)];
    }

    /*!
     * \brief Return the SAMRAI variable represented by the fluid variable name.
     */
    inline const std::string& getFLName(const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>& fl_var)
    {
        auto fl_it = std::find(d_fl_vars.begin(), d_fl_vars.end(), fl_var);
        TBOX_ASSERT(fl_it != d_fl_vars.end());
        return d_fl_names[std::distance(d_fl_vars.begin(), fl_it)];
    }

    /*!
     * \brief Get the pair of coupling lists for the given surface quantity. The first vector contains the surface
     * system names and the second vector contains the fluid system names.
     */
    inline void getSFCouplingLists(const std::string& sf_name,
                                   std::vector<std::string>& sf_names,
                                   std::vector<std::string>& fl_names)
    {
        sf_names = d_sf_sf_map[sf_name];
        fl_names = d_sf_fl_map[sf_name];
    }

    /*!
     * \brief Get the pair of coupling lists for the given fluid quantity. The first vector contains the surface system
     * names and the second vector contains the fluid system names.
     */
    inline void getFLCouplingLists(const std::string& fl_name,
                                   std::vector<std::string>& sf_names,
                                   std::vector<std::string>& fl_names)
    {
        sf_names = d_fl_sf_map[fl_name];
        fl_names = d_fl_fl_map[fl_name];
    }

    /*!
     * \brief Return the reaction function and function context pair for the given surface quantity.
     */
    inline const ReactionFcnCtx& getSFReactionFcnCtxPair(const std::string& sf_name)
    {
        return d_sf_reaction_fcn_ctx_map[sf_name];
    }

    /*!
     * \brief Return the boundary condition functions for the given fluid concentration. Note that the first element is
     * the a function and the second element is the g function.
     *
     * @{
     */
    inline const BdryConds&
    getFLBdryConditionFcns(const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>& fl_var)
    {
        return getFLBdryConditionFcns(fl_var->getName());
    }
    inline const BdryConds& getFLBdryConditionFcns(const std::string& fl_name)
    {
        return d_fl_a_g_fcn_map[fl_name];
    }
    /*!
     * @}
     */

private:
    std::string d_object_name;
    libMesh::Mesh* d_mesh = nullptr;
    IBTK::FEDataManager* d_fe_data_manager = nullptr;

    std::shared_ptr<CutCellMeshMapping> d_cut_cell_mesh_mapping = nullptr;

    std::vector<std::string> d_sf_names;
    std::vector<std::string> d_fl_names;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>> d_fl_vars;

    std::map<std::string, std::vector<std::string>> d_sf_fl_map;
    std::map<std::string, std::vector<std::string>> d_sf_sf_map;
    std::map<std::string, std::vector<std::string>> d_fl_fl_map;
    std::map<std::string, std::vector<std::string>> d_fl_sf_map;
    std::map<std::string, ReactionFcnCtx> d_sf_reaction_fcn_ctx_map;
    std::map<std::string, BdryConds> d_fl_a_g_fcn_map;

    std::string d_J_sys_name;

    RBFReconstructCache d_rbf_reconstruct;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    int d_stencil_size = -1;
    bool d_update_weights = false;

    bool d_fe_eqs_initialized = false;

    int d_vol_idx = IBTK::invalid_index, d_ls_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_scr_var;
    int d_scr_idx = IBTK::invalid_index;
};

} // namespace LS
#endif
