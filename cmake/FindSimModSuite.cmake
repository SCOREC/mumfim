# - Try to find Simmetrix SimModSuite
# Once done this will define
#  SIMMODSUITE_FOUND - System has SimModSuite
#  SIMMODSUITE_INCLUDE_DIRS - The SimModSuite include directories
#  SIMMODSUITE_LIBS - The libraries needed to use SimModSuite
#  SIMMODSUITE_DEFINITIONS - Compiler switches required for using SimModSuite
#

find_library(simpartitioned SimPartitionedMesh-mpi)
find_library(simmeshing SimMeshing)
find_library(simmodel SimModel)
find_library(simmeshtools SimMeshTools)
find_library(simwrapper SimPartitionWrapper-${SIM_MPI})
find_library(simfield SimField)
find_library(simparasolid SimParasolid270)
set(SIMMODSUITE_LIBS "${simmeshing};${simfield};${simmeshtools};${simpartitioned};${simwrapper};${simparasolid};${simmodel}")

find_path(SIMMODSUITE_INCLUDE_DIR 
  NAMES SimUtil.h SimError.h SimModel.h 
  PATHS ${SIM_INCLUDE_DIR})
if(NOT EXISTS "${SIMMODSUITE_INCLUDE_DIR}")
  message(STATUS "simmetrix include dir not found")
endif()

set(SIMMODSUITE_INCLUDE_DIRS ${SIMMODSUITE_INCLUDE_DIR})

string(REGEX REPLACE 
  "/include$" "" 
  SIMMODSUITE_INSTALL_DIR
  "${SIMMODSUITE_INCLUDE_DIR}")

string(REGEX MATCH 
  "[0-9]+.[0-9]+-[0-9]+"
  SIM_VERSION
  "${SIMMODSUITE_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SIMMODSUITE  DEFAULT_MSG
                                  SIMMODSUITE_LIBS SIMMODSUITE_INCLUDE_DIRS)

mark_as_advanced(SIMMODSUITE_INCLUDE_DIR SIMMODSUITE_LIBS)

set(SIM_LINK_LIBS "")
foreach(lib ${SIM_LIB_NAMES})
  set(SIM_LINK_LIBS "${SIM_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${SIMMODSUITE_INSTALL_DIR}")
set(includedir "${SIMMODSUITE_INCLUDE_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libSimModSuite.pc.in"
  "${CMAKE_BINARY_DIR}/libSimModSuite.pc"
  @ONLY)

#is this OK for a package file????
if(NOT BUILD_IN_TRILINOS)
  INSTALL(FILES "${CMAKE_BINARY_DIR}/libSimModSuite.pc" DESTINATION lib/pkgconfig)
endif()
