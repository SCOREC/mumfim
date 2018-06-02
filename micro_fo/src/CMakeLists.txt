set(micro_fo_public_headers
  bioFiberNetwork.h
  bioFiberRVEAnalysis.h
  bioFiberRVEAnalysis_impl.h
  bioMicroFOMultiscale.h
  bioMultiscaleRVEAnalysis.h
  bioMultiscaleRVE.h
  bioRVE.h
  bioRVEVolumeTerms.h
  bioTrussIntegrator.h
  bioUtil.h
  bioUtil_impl.h
  lasCSR.h
  lasSparskit.h
  lasSparskitExterns.h)

set(micro_fo_src_files
  bioRVE.cc
  bioFiberNetwork.cc
  bioFiberNetworkIO.cc
  bioFiberRVEAnalysis.cc
  bioAlignFiberNetwork.cc
  bioMicroFOMultiscale.cc
  bioMultiscaleRVEAnalysis.cc
  bioMultiscaleRVE.cc
  bioUtil.cc
  lasCSR.cc
  lasCSRBuilder.cc
  lasSparskit.cc
  lasSparskitExterns.cc)

find_package(sparskit REQUIRED)
find_package(amsi COMPONENTS util analysis multiscale REQUIRED)

add_library(micro_fo ${micro_fo_src_files} ${micro_fo_public_headers})
set_property(TARGET micro_fo PROPERTY VERSION ${BIOTISSUE_VERSION})

target_include_directories(micro_fo PUBLIC
                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                    $<INSTALL_INTERFACE:include/biotissue>)

target_link_libraries(micro_fo PUBLIC
                      ${SPARSKIT_LIBRARY}
                      amsi::util amsi::analysis amsi::multiscale)

install(TARGETS micro_fo
        EXPORT biotissueTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include/biotissue)

install(FILES ${micro_fo_header_files} DESTINATION include/biotissue)

