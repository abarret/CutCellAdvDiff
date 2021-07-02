ADD_LIBRARY(CCAD2d)
ADD_LIBRARY(CCAD3d)

SET(FORTRAN_SRC_BASE
	fortran/advection_routines)

FOREACH(_d ${CCAD_DIMENSIONS})
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    CCAD_PROCESS_M4(NAME ${_fortran_root_name}${_d}d.f)
  ENDFOREACH()

  SET("FORTRAN_GENERATED_SRC${_d}D" "")
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    LIST(APPEND "FORTRAN_GENERATED_SRC${_d}D" "${CMAKE_CURRENT_BINARY_DIR}/${_fortran_root_name}${_d}d.f")
  ENDFOREACH()
ENDFOREACH()


SET(CXX_SRC
    AdvectiveReconstructionOperator.cpp
    CutCellMeshMapping.cpp
    CutCellVolumeMeshMapping.cpp
    FEMeshPartitioner.cpp
    IntegrateFunction.cpp
    LSAdvDiffIntegrator.cpp
    LSCartGridFunction.cpp
    LSCutCellBoundaryConditions.cpp
    LSCutCellLaplaceOperator.cpp
    LSFindCellVolume.cpp
    LSFromLevelSet.cpp
    LSFromMesh.cpp
    ls_functions.cpp
    MLSReconstructCache.cpp
    RBFReconstructCache.cpp
    RBFReconstructions.cpp
    ReconstructCache.cpp
    reconstructions.cpp
    SBAdvDiffIntegrator.cpp
    SBBoundaryConditions.cpp
    SBIntegrator.cpp
    SBSurfaceFluidCouplingManager.cpp
    VolumeBoundaryMeshMapping.cpp
    ZSplineReconstructions.cpp
    )

TARGET_SOURCES(CCAD2d PRIVATE ${FORTRAN_GENERATED_SRC2D} ${CXX_SRC})
TARGET_SOURCES(CCAD3d PRIVATE ${FORTRAN_GENERATED_SRC3D} ${CXX_SRC})

TARGET_INCLUDE_DIRECTORIES(CCAD2d PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include/>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>)
TARGET_INCLUDE_DIRECTORIES(CCAD3d PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include/>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>)

TARGET_LINK_LIBRARIES(CCAD2d PUBLIC CCADHeaders)
TARGET_LINK_LIBRARIES(CCAD3d PUBLIC CCADHeaders)

CCADLIB_SETUP_TARGET_LIBRARY(CCAD2d)
CCADLIB_SETUP_TARGET_LIBRARY(CCAD3d)

INSTALL(TARGETS CCAD2d EXPORT CCADTargets COMPONENT library)
INSTALL(TARGETS CCAD3d EXPORT CCADTargets COMPONENT library)
