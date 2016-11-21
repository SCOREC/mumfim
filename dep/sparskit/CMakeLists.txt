cmake_minimum_required (VERSION 2.6)
project (sparskit_scorec)
set(PACKAGE_VERSION "1.0")

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

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

endif()
# mpi }

# Add local directory
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src)

# Include LAS headers
include_directories(${LAS_INCLUDE_DIR})

set(HEADER_FILES 
  src/CSR_Matrix.h
  src/CSR_Vector.h
  src/listman.h)
#  src/SparskitLinearSystemSolver.h )

set(SOURCE_FILES 
  src/CSR_Matrix.cc
  src/listman.cc
  src/CSR_SPARSKIT_Solver.cc )
#  src/SparskitLinearSystemSolver.cc
#  src/CSR_SPARSKIT_Solver.cc )

# Add the library to create
add_library(SPARSKIT_SCOREC ${SOURCE_FILES} ${HEADER_FILES})

# Link with mpi libraries
target_link_libraries(SPARSKIT_SCOREC ${MPI_LIBRARIES})

# For install, put the created AMSI library in installpath/lib
install (TARGETS SPARSKIT_SCOREC DESTINATION lib)

# For install, put headers in installpath/include
install (FILES ${HEADER_FILES} DESTINATION  include)

# configure the modulefile
string(TOUPPER ${PROJECT_NAME} PACKAGE_NAME)
set(MODULEFILE ${PROJECT_NAME}_${PACKAGE_VERSION})
configure_file(moduletemplate ${MODULEFILE} @ONLY)

# install the modulefile 
install(FILES build/${MODULEFILE} DESTINATION lib/module)



