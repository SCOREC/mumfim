set(micro_fo_src_files
  bioBatchedFiberRVEAnalysisExplicit.cc
  bioExplicitOutputWriter.cc
  bioFiberNetwork.cc
  bioFiberNetworkIO.cc
  bioFiberReactions.cc
  bioFiberRVEAnalysis.cc
  bioFiberRVEAnalysisCreation.cc
  bioFiberRVEAnalysisExplicit.cc
  bioFiberRVEAnalysisStaticImplicit.cc
  bioMassIntegrator.cc
  bioMicroFOParams.cc
  bioMultiscaleCoupling.cc
  bioMultiscaleMicroFOParams.cc
  bioMultiscaleRVEAnalysis.cc
  bioNeoHookeanRVEAnalysis.cc
  bioRVE.cc
  bioRVEAnalysis.cc
  bioTrussIntegrator.cc
  bioUtility.cc
  bioFiberNetworkLibrary.cc
  bioApfPointers.cc)

find_package(las HINTS "$ENV{DEVROOT}/install/las/lib/cmake/" REQUIRED) #COMPONENTS sparskit core)
find_package(las_core HINTS "$ENV{DEVROOT}/install/las/lib/cmake" REQUIRED)
find_package(amsi COMPONENTS util analysis multiscale REQUIRED)
find_package(yaml-cpp REQUIRED)

find_package(Kokkos CONFIG REQUIRED)

option(ENABLE_WRITE_MICRO_PER_ITER "Write microscale to vtk on every iteration" OFF)
if(ENABLE_WRITE_MICRO_PER_ITER)
  set(VERBOSITY_FLAGS "${VERBOSITY_FLAGS} -DWRITE_MICRO_PER_ITER")
endif()

option(ENABLE_WRITE_MICRO_PER_STEP "Write microscale to vtk on every step" OFF)
if(ENABLE_WRITE_MICRO_PER_STEP)
  set(VERBOSITY_FLAGS "${VERBOSITY_FLAGS} -DWRITE_MICRO_PER_STEP")
endif()

add_library(micro_fo ${micro_fo_src_files})
set_property(TARGET micro_fo PROPERTY VERSION ${BIOTISSUE_VERSION})
if(VERBOSITY_FLAGS)
    set_target_properties(micro_fo PROPERTIES COMPILE_FLAGS ${VERBOSITY_FLAGS})
endif()

target_include_directories(micro_fo PUBLIC
                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/batched_analysis>
                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                                    $<INSTALL_INTERFACE:include/biotissue>
                                    PRIVATE ${yaml-cpp_INCLUDE_DIRS} )

target_link_libraries(micro_fo PUBLIC las amsi::util amsi::analysis amsi::multiscale 
                      las_core Kokkos::kokkos
                      PRIVATE ${yaml-cpp_LIBRARIES})
if(yaml-cpp_FOUND)
  target_compile_definitions(micro_fo PUBLIC -DHAS_YAML_CPP)
endif(yaml-cpp_FOUND)

set(ENABLE_VERBOSITY HIGH CACHE STRING "Enable Debug/Verbosity Level")
set_property(CACHE ENABLE_VERBOSITY PROPERTY STRINGS HIGH MED LOW OFF)
if(ENABLE_VERBOSITY)
    set(BIO_VERBOSITY_DEFINITIONS "-DBIO_VERBOSE_1")
    if(ENABLE_VERBOSITY STREQUAL "LOW")
    set(BIO_VERBOSITY_DEFINITIONS "-DBIO_VERBOSE_1")
  elseif(ENABLE_VERBOSITY STREQUAL "MED")
    set(BIO_VERBOSITY_DEFINITIONS "-DBIO_VERBOSE_1 -DBIO_VERBOSE_2")
  elseif(ENABLE_VERBOSITY STREQUAL "HIGH")
    set(BIO_VERBOSITY_DEFINITIONS "-DBIO_VERBOSE_1 -DBIO_VERBOSE_2 -DBIO_VERBOSE_3")
  elseif(ENABLE_VERBOSITY STREQUAL "OFF")
    set(BIO_VERBOSITY_DEFINITIONS "")
  else()
    message(FATAL_ERROR "Something went wrong with VERBOSITY setting: ${ENABLE_VERBOSITY}\n Valid options are LOW,MED,HIGH,OFF")
  endif()
endif(ENABLE_VERBOSITY)
target_compile_definitions(micro_fo PUBLIC ${BIO_VERBOSITY_DEFINITIONS})

set(ENABLE_VERBOSITY HIGH CACHE STRING "Enable Debug/Verbosity Level")
set_property(CACHE ENABLE_VERBOSITY PROPERTY STRINGS HIGH MED LOW OFF)
if(ENABLE_VERBOSITY)
    set(BIO_VERBOSITY_DEFINITIONS "-DBIO_VERBOSE_1")
    if(ENABLE_VERBOSITY STREQUAL "LOW")
    set(BIO_VERBOSITY_DEFINITIONS "-DBIO_VERBOSE_1")
  elseif(ENABLE_VERBOSITY STREQUAL "MED")
    set(BIO_VERBOSITY_DEFINITIONS "-DBIO_VERBOSE_1 -DBIO_VERBOSE_2")
  elseif(ENABLE_VERBOSITY STREQUAL "HIGH")
    set(BIO_VERBOSITY_DEFINITIONS "-DBIO_VERBOSE_1 -DBIO_VERBOSE_2 -DBIO_VERBOSE_3")
  elseif(ENABLE_VERBOSITY STREQUAL "OFF")
    set(BIO_VERBOSITY_DEFINITIONS "")
  else()
    message(FATAL_ERROR "Something went wrong with VERBOSITY setting: ${ENABLE_VERBOSITY}\n Valid options are LOW,MED,HIGH,OFF")
  endif()
endif(ENABLE_VERBOSITY)
target_compile_definitions(micro_fo PUBLIC ${BIO_VERBOSITY_DEFINITIONS})

set(MICRO_BACKEND sparskit CACHE STRING "sparskit petsc")
if(MICRO_BACKEND STREQUAL "sparskit")
  target_compile_definitions(micro_fo PUBLIC -DMICRO_BACKEND=sparskit -DMICRO_USING_SPARSKIT)
elseif(MICRO_BACKEND STREQUAL "petsc")
  target_compile_definitions(micro_fo PUBLIC -DMICRO_BACKEND=petsc -DMICRO_USING_PETSC)
else()
  message(FATAL_ERROR "You have incorrectly set MICRO_BACKEND=${MICRO_BACKEND}.\n Balid options are sparskit, petsc")
endif()

#target_include_directories(micro_fo SYSTEM PRIVATE $ENV{INCLUDE})

install(TARGETS micro_fo
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include/biotissue)
