ADD_LIBRARY(LSLib2d)
ADD_LIBRARY(LSLib3d)

SET(FORTRAN_SRC_BASE
	fortran/advection_routines)

FOREACH(_d ${IBAMR_DIMENSIONS})
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    IBAMR_PROCESS_M4(NAME ${_fortran_root_name}${_d}d.f)
  ENDFOREACH()

  SET("FORTRAN_GENERATED_SRC${_d}D" "")
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    LIST(APPEND "FORTRAN_GENERATED_SRC${_d}D" "${CMAKE_CURRENT_BINARY_DIR}/${_fortran_root_name}${_d}d.f")
  ENDFOREACH()
ENDFOREACH()


SET(CXX_SRC
    IntegrateFunction.cpp
    LSCartGridFunction.cpp
    LSCutCellBoundaryConditions.cpp
    LSCutCellLaplaceOperator.cpp
    LSFindCellVolume.cpp
    LSFromLevelSet.cpp
    LSFromMesh.cpp
    QInitial.cpp
    RBFReconstructCache.cpp
    SBBoundaryConditions.cpp
    SBIntegrator.cpp
    SemiLagrangianAdvIntegrator.cpp
    SetLSValue.cpp
    )

TARGET_SOURCES(LSLib2d PRIVATE ${FORTRAN_GENERATED_SRC2D} ${CXX_SRC})
TARGET_SOURCES(LSLib3d PRIVATE ${FORTRAN_GENERATED_SRC3D} ${CXX_SRC})

TARGET_LINK_LIBRARIES(LSLib2d PUBLIC LSHeaders)
TARGET_LINK_LIBRARIES(LSLib3d PUBLIC LSHeaders)

LS_SETUP_TARGET_LIBRARY(LSLib2d)

INSTALL(TARGETS LSLib2d EXPORT LSTargets COMPONENT library)
INSTALL(TARGETS LSLib3d EXPORT LSTargets COMPONENT library)
