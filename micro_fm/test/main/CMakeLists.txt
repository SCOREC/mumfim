cmake_minimum_required (VERSION 2.6)
project (MICRO_FM)

message("Building on " ${HOSTNAME} " system...")

if(NOT ${MANUAL_SET_MPI})
  message("Using find_package(MPI)")
  find_package(MPI REQUIRED)
else(NOT ${MANUAL_SET_MPI})
  message("Not using find_package(MPI)")
endif(NOT ${MANUAL_SET_MPI})

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
SET(CMAKE_C_COMPILER ${MPI_COMPILER})
SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})

# Add mpi link flags
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS})

# Debug or release
#if($ENV{BIOTISSUE_OPTIMIZE})
#  set(CMAKE_BUILD_TYPE Release)
#else()
  set(CMAKE_BUILD_TYPE Debug)
#endif()

set(SOURCE_FILES main.cc)

add_executable(MICRO_FM main.cc)

link_directories(${PETSC_PATH}/${PETSC_ARCH}/lib
	/usr/X11/R6/lib64
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/psKrnl)

# Set libraries
set(MICRO_FM_LIBRARIES 
	${BIOTISSUE_LIBRARY_PATH}/libMICRO_FM.a
	${BIOTISSUE_LIBRARY_PATH}/libFEA.a
	${BIOTISSUE_LIBRARY_PATH}/libLAS.a
	${BIOTISSUE_LIBRARY_PATH}/libapf.a
	${BIOTISSUE_LIBRARY_PATH}/libapf_sim.a
	${BIOTISSUE_LIBRARY_PATH}/libpcu.a
	${BIOTISSUE_LIBRARY_PATH}/libAMSI.a
	${PETSC_PATH}/${PETSC_ARCH}/lib/libpetsc.a
	${PETSC_PATH}/${PETSC_ARCH}/lib/libsuperlu_dist_2.5.a
	${PETSC_PATH}/${PETSC_ARCH}/lib/libparmetis.a
	${PETSC_PATH}/${PETSC_ARCH}/lib/libmetis.a
	${PETSC_PATH}/${PETSC_ARCH}/lib/libflapack.a
	${PETSC_PATH}/${PETSC_ARCH}/lib/libfblas.a 
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/libSimMeshing.a
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/libSimMeshTools.a
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/libSimPartitionedMesh-mpi.a
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/libSimModel.a
        )

# Set system specific libraries

# BlueGeneQ
if(${HOSTNAME} STREQUAL bgq)
  set(MICRO_FM_LIBRARIES ${MICRO_FM_LIBRARIES}
        ${BIOTISSUE_LIBRARY_PATH}/libpcu.a
        ${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/libSimPartitionWrapper-bgmpi.a
	/opt/ibmcmp/vacpp/bg/12.1/lib64/libibmc++.a
	/opt/ibmcmp/vacpp/bg/12.1/lib64/libxl.a
	/bgsys/drivers/toolchain/V1R2M0-efix23/gnu-linux/powerpc64-bgq-linux/lib/libdl.a
	/usr/lib/gcc/ppc64-redhat-linux/4.4.4/libgfortran.a
	/bgsys/drivers/toolchain/V1R2M0-efix23/gnu-linux/powerpc64-bgq-linux/lib/libpthread.a
	/gpfs/bgdata/bgsys/drivers/toolchain/V1R2M0-efix23/gnu-linux/powerpc64-bgq-linux/lib/libm.a
        )
endif(${HOSTNAME} STREQUAL bgq)

# SCOREC system
if(${HOSTNAME} STREQUAL scorec)
  set(MICRO_FM_LIBRARIES ${MICRO_FM_LIBRARIES}
        ${BIOTISSUE_LIBRARY_PATH}/libpcu.a
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/libSimPartitionWrapper-openmpi.a
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/libSimParasolid260.a
	${MESHSIM_INCLUDE_PATH}/lib/${MESHSIM_LIB_SUBDIR}/psKrnl/libpskernel.so
	-lX11
	-lgfortran
	-lpthread
	-lm
        )
endif(${HOSTNAME} STREQUAL scorec)

target_link_libraries(MICRO_FM ${MICRO_FM_LIBRARIES})

set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS})

# Set includes
include_directories(${PETSC_PATH}/include)
include_directories(${PETSC_PATH}/${PETSC_ARCH}/include)
include_directories(${BOOST_INCLUDE_PATH})
include_directories(${AMSI_INCLUDE_PATH})
include_directories(${MESHSIM_INCLUDE_PATH}/include)
include_directories(${FEA_INCLUDE_PATH})


