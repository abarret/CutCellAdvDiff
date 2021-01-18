CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

PROJECT(LS
  DESCRIPTION "Advection-Diffusion solver with reaction boundary conditions"
  LANGUAGES C CXX Fortran)

SET(LS_DIMENSIONS "2" "3")
FIND_PACKAGE(IBAMR REQUIRED HINTS ${IBAMR_DIR} $ENV{IBAMR_DIR})
  
FUNCTION(LS_SETUP_TARGET_LIBRARY target_library)
	MESSAGE(STATUS "setting up target ${target_library}")

	SET(_2_location = -1)
	SET(_3_location = -1)
	STRING(FIND "${target_library}" "2" _2_location)
	STRING(FIND "${target_library}" "3" _3_location)
	SET(_d "")
	IF(NOT ${_2_location} STREQUAL "-1")
		SET(_d "2")
	ELSE()
		SET(_d "3")
	ENDIF()

	MESSAGE(STATUS "Adding flag -DNDIM=${_d} to target ${target_library}")
	TARGET_COMPILE_OPTIONS(${target_library} PUBLIC -DNDIM=${_d})
	TARGET_LINK_LIBRARIES(${target_library} PUBLIC IBAMR::IBAMR${_d}d)
ENDFUNCTION()

ADD_LIBRARY(LSHeaders INTERFACE)
TARGET_INCLUDE_DIRECTORIES(LSHeaders INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>)

INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/LS DESTINATION include)

INSTALL(TARGETS LSHeaders EXPORT LSTargets)

ADD_SUBDIRECTORY(src)
#ADD_SUBDIRECTORY(examples)
