CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

PROJECT(CCAD LANGUAGES C CXX Fortran)
SET(CCAD_VERSION "0.1.0")
SET(CCAD_DIMENSIONS
    "2"
    "3")

OPTION(BUILD_SHARED_LIBS "Whether or not to build shared libraries." ON)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

MESSAGE(STATUS "This is CMAKE ${CMAKE_VERSION}")
MESSAGE(STATUS "")
MESSAGE(STATUS "This is CCAD Version ${CCAD_VERSION}")
INCLUDE(GNUInstallDirs)
INCLUDE(CMakePackageConfigHelpers)
LIST(APPEND CMAKE_MODULE_PARTH ${CMAKE_SOURCE_DIR}/cmake/modules)

MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up IBAMR")
FIND_PACKAGE(IBAMR REQUIRED)

INCLUDE(CMakeParseArguments)
#
# Macro to process m4 files to generate Fortran.
#
# This macro includes some ad - hoc logic to attempt to detect M4 dependencies
# (i.e., include() statements in the provided file) so that the build system is
# aware of them.However, since this is done with regular expressions it is
# fragile and may not always pick up file dependencies.For example, writing
#
# input(
# dnl skip this line
# / path / to / file.m4
# )
#
# is perfectly valid M4 but beyond the capacity of this macro to parse.In such
# cases we simply do not add the dependency and print a warning - compilation
# from scratch will still succeed but the build system will not be able to track
# the dependency.
#
# Required arguments:
#
# - NAME : Relative path(from the directory in which the macro is called) to the
# relevant Fortran file.For example : if we want to preprocess
# foo / bar / baz.f.m4 then we should pass foo / bar / baz.f as the NAME argument.
MACRO(CCAD_PROCESS_M4)
SET(options)
SET(singleValueArgs NAME)

CMAKE_PARSE_ARGUMENTS(arg
                      "${options}"
                      "${singleValueArgs}"
                      "${multiValueArgs}" ${ARGN})
GET_FILENAME_COMPONENT(_in_directory "${CMAKE_CURRENT_SOURCE_DIR}/${arg_NAME}.m4" DIRECTORY)
GET_FILENAME_COMPONENT(_out_directory "${CMAKE_CURRENT_BINARY_DIR}/${arg_NAME}.f" DIRECTORY)
FILE(MAKE_DIRECTORY ${_out_directory})

SET(_samrai_fortdir ${SAMRAI_ROOT}/include)
SET(_current_srcdir ${_in_directory})
SET(_top_srcdir ${CMAKE_SOURCE_DIR})
SET(_input "${CMAKE_CURRENT_SOURCE_DIR}/${arg_NAME}.m4")
SET(_output "${CMAKE_CURRENT_BINARY_DIR}/${arg_NAME}")
SET(_args
    "-DSAMRAI_FORTDIR=${_samrai_fortdir}"
    "-DCURRENT_SRCDIR=${_current_srcdir}"
    "-DTOP_SRCDIR=${_top_srcdir}")

FILE(STRINGS ${_input} _parsed_includes REGEX "include.*DIR")
SET(_expanded_files)
FOREACH(_parsed_include ${_parsed_includes})
STRING(REGEX REPLACE
       "include\\((.*)\\)(dnl)?$"
       "\\1" _tmp ${ _parsed_include })
STRING(REGEX REPLACE
       "SAMRAI_FORTDIR"
       "${_samrai_fortdir}" _tmp ${_tmp})
STRING(REGEX REPLACE
       "CURRENT_SRCDIR"
       "${_current_srcdir}" _tmp ${_tmp})
STRING(REGEX REPLACE
       "TOP_SRCDIR"
       "${_top_srcdir}" _tmp ${_tmp})
LIST(APPEND _expanded_files ${_tmp})
ENDFOREACH()

SET(_dependencies)
LIST(APPEND _dependencies ${_input})
FOREACH(_expanded_file ${_expanded_files})
IF(EXISTS ${_expanded_file})
LIST(APPEND _dependencies ${_expanded_file})
ELSE()
      MESSAGE(WARNING
"Unable to locate file\n    ${_expanded_file}\n which is a dependency of
${_input_file}. Since dependency tracking of M4 files is done manually this is
not a fatal problem but it may result in compilation failures later in the build process.")
    ENDIF()
  ENDFOREACH()

  ADD_CUSTOM_COMMAND(
    OUTPUT ${_output}
    DEPENDS ${_dependencies}
    COMMAND m4 ${_args} ${_input} > ${_output}
    VERBATIM)

  SET_SOURCE_FILES_PROPERTIES(${_output} PROPERTIES GENERATED true)
  SET_SOURCE_FILES_PROPERTIES(${_output} PROPERTIES LANGUAGE Fortran)
ENDMACRO()

# Macro to setup an CCAD target library with all common features(i.e.,
# everything but the source files)
FUNCTION(CCADLIB_SETUP_TARGET_LIBRARY target_library)
  MESSAGE(STATUS "setting up target ${target_library}")

  SET_TARGET_PROPERTIES(${target_library} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

# Figure out if we are a 2D or a 3D library:
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

#
# CCADLib is compiled in 2D and 3D, but the headers are only installed
# once - hence set up the headers as separate targets
#
ADD_LIBRARY(CCADHeaders INTERFACE)
TARGET_INCLUDE_DIRECTORIES(
  CCADHeaders
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>)

INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/CCAD DESTINATION include)

INSTALL(TARGETS CCADHeaders EXPORT CCADTargets)

#
# Set up specific targets for executables and libraries that want to link to
# CCAD.As noted above, we have to do this since SAMRAI might not be compiled
# with - fPIC : i.e., if we added SAMRAI as a link dependency to libIBTK and
# libIBTK as a link dependency to libFoo, then the linker would encounter errors
# as a result.With executables we can use the full link interface.
#
SET(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ccad/)

# Set up actual files containing the export target information:
INSTALL(EXPORT "CCADTargets" FILE "CCADTargets.cmake"
	NAMESPACE CCAD::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ccad/)

# Do the same for in - tree:
EXPORT(EXPORT "CCADTargets"
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/CCADTargets.cmake"
  NAMESPACE CCAD::)

WRITE_BASIC_PACKAGE_VERSION_FILE(
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/CCADConfigVersion.cmake"
  VERSION ${CCAD_VERSION}
  COMPATIBILITY AnyNewerVersion)

CONFIGURE_PACKAGE_CONFIG_FILE(
  ${CMAKE_SOURCE_DIR}/cmake/CCADConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CCADConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/CCADConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CCADConfigVersion.cmake
  DESTINATION ${INSTALL_CONFIGDIR})

# proceed to compilation units:
ADD_SUBDIRECTORY(src)

ADD_SUBDIRECTORY(examples)
