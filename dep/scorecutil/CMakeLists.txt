cmake_minimum_required (VERSION 2.6)
project (SCORECUTIL)
set (PACKAGE_VERSION "1.0")

# mpi {
if(${MANUAL_SET_MPI})

# Output all mpi associated variables
message("MPI_FOUND = " ${MPI_FOUND})
message("MPI_INCLUDE_PATH = " ${MPI_INCLUDE_PATH})
message("MPI_COMPILE_FLAGS = " ${MPI_COMPILE_FLAGS})
message("MPI_COMPILER = " ${MPI_COMPILER})
message("MPI_LINK_FLAGS = " ${MPI_LINK_FLAGS})
message("MPI_LIBRARIES = " ${MPI_LIBRARIES})
message("MPIEXEC = " ${MPIEXEC})
message("MPIEXEC_NUMPROC_FLAG = " ${MPIEXEC_NUMPROC_FLAG})

# Set compilers to mpi
set(CMAKE_C_COMPILER ${MPI_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_COMPILER})

# Add mpi link flags
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS})

else()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpicxx)
else()
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

endif()

# Add simmetrix headers
include_directories(${SIM_INCLUDE_DIR})

# Add boost headers
include_directories(${BOOST_INCLUDE_DIR})

# Add local include directory
include_directories(include)
set(HEADER_FILES 
  include/affine_space.h
  include/BezierMappingBuilder.h
  include/BezierMapping.h
  include/GaussIntegrator.h
  include/GaussLegendreSimplex.h
  include/Integrator.h
  include/GaussQuadrature.h
  include/IntPt.h
  include/LagrangeMappingBuilder.h
  include/LagrangeMapping.h
  include/MappingBuilder.h
  include/Mapping.h
  include/mPoint.h
  include/MSMapping.h
  include/mTensor2.h
  include/mTensor2Symmetric.h
  include/mTensor4.h
  include/mVector.h
  include/OctreeCreate2.h
  include/Octree.h
  include/o_internals.h
  include/scorec_function_objects.h
  include/scorecMD_memory.h
  include/scorecSC_mem.h
  include/scorecSListBase.h
  include/scorecSSList.h
  include/CountTime.h
  include/VolumeBuckets.h)

set(SOURCE_FILES 
  calculus/BezierMappingBuilder.cc
  calculus/BezierMapping.cc
  calculus/GaussLegendre1D.cc
  calculus/GaussLegendreSimplex.cc
  calculus/GaussQuadrature.cc
  calculus/GaussQuadratureHex.cc
  calculus/GaussQuadratureQuad.cc
  calculus/GaussQuadratureTet.cc
  calculus/GaussQuadratureTri.cc
  calculus/GaussQuadraturePrism.cc
  calculus/GaussQuadraturePyramid.cc
  calculus/Integrator.cc
  calculus/LagrangeMappingBuilder.cc
  calculus/LagrangeMapping.cc
  calculus/Mapping.cc
  octree/OctreeSrc.cc
  octree/o_internals.cc
  util/scorecMD_memory.cc
  util/scorecSC_mem.cc
  util/scorecSListBase.cc)

# project-specific
SET_SOURCE_FILES_PROPERTIES(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "-DSIM")

# Add the library to create
add_library(SCORECUtil ${SOURCE_FILES} ${HEADER_FILES})

# Link with mpi libraries
target_link_libraries(SCORECUtil ${MPI_LIBRARIES})

# For install, put the created AMSI library in installpath/lib
install (TARGETS SCORECUtil DESTINATION lib)

# For install, put headers in installpath/include
install (FILES ${HEADER_FILES} DESTINATION  include)

# modulefile configure and install
string(TOUPPER ${PROJECT_NAME} PACKAGE_NAME)
set(MODULEFILE ${PROJECT_NAME}_${PACKAGE_VERSION})
configure_file(moduletemplate ${MODULEFILE} @ONLY)

# For install, put the created modulefile in the install path
install(FILES build/${MODULEFILE} DESTINATION lib/module)



