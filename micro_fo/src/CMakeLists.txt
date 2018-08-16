set(micro_fo_public_headers
  bioFiberNetwork.h
  bioFiberRVEAnalysis.h
  bioFiberRVEAnalysis_impl.h
  bioMultiscaleCoupling.h
  bioMultiscaleMicroFOParams.h
  bioMultiscaleRVEAnalysis.h
  bioMultiscaleRVE.h
  bioRVE.h
  bioRVEVolumeTerms.h
  bioTrussIntegrator.h
  bioUtil.h
  bioUtil_impl.h
  bioVerbosity.h)

set(micro_fo_src_files
  bioRVE.cc
  bioFiberNetwork.cc
  bioFiberNetworkIO.cc
  bioFiberRVEAnalysis.cc
  bioAlignFiberNetwork.cc
  bioMultiscaleCoupling.cc
  bioMultiscaleMicroFOParams.cc
  bioMultiscaleRVEAnalysis.cc
  bioMultiscaleRVE.cc
  bioUtil.cc)

#find_package(sparskit REQUIRED)
find_package(las HINTS "$ENV{DEVROOT}/install/las/lib/cmake/" REQUIRED) #COMPONENTS sparskit core)
find_package(las_core HINTS "$ENV{DEVROOT}/install/las/lib/cmake" REQUIRED)
find_package(amsi COMPONENTS util analysis multiscale REQUIRED)

if(WRITE_MICRO_ITER)
  set(VERBOSITY_FLAGS "${VERBOSITY_FLAGS} -DWRITE_MICRO_PER_ITER")
endif()

if(WRITE_MICRO_STEP)
  set(VERBOSITY_FLAGS "${VERBOSITY_FLAGS} -DWRITE_MICRO_PER_STEP")
endif()

add_library(micro_fo ${micro_fo_src_files} ${micro_fo_public_headers})
set_property(TARGET micro_fo PROPERTY VERSION ${BIOTISSUE_VERSION})
if(VERBOSITY_FLAGS)
    set_target_properties(micro_fo PROPERTIES COMPILE_FLAGS ${VERBOSITY_FLAGS})
endif()

target_include_directories(micro_fo PUBLIC
                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                    $<INSTALL_INTERFACE:include/biotissue>)

target_link_libraries(micro_fo PUBLIC las amsi::util amsi::analysis amsi::multiscale las_core)

install(TARGETS micro_fo
        EXPORT biotissueTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include/biotissue)

install(FILES ${micro_fo_header_files} DESTINATION include/biotissue)

