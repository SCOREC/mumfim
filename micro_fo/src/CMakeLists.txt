set(micro_fo_public_headers
  bioRVE.h
  bioRVE2.h
  bioFiberNetwork.h
  bioFiberRVEAnalysis.h
  bioFiberRVEAnalysis_impl.h
  bioRVEUtil.h
  bioMicroFOMultiscale.h
  bioMultiscaleRVEAnalysis.h
  bioMultiscaleRVE.h
  bioTrussIntegrator.h
  bioUtil.h
  bioUtil_impl.h
  lasCSR.h
  lasSparskit.h
  lasSparskitExterns.h)

set(micro_fo_src_files
  bioRVE.cc
  bioRVE2.cc
  bioFiberNetwork.cc
  bioFiberNetworkIO2.cc
  bioFiberRVEAnalysis.cc
  bioAlignFiberNetwork.cc
  bioSolverPsi.cc
  bioSolverGMRES.cc
  bioSolverBeam.cc
  bioMacroElement.cc
  bioMacroTerms.cc
  bioMultiscaleRVEAnalysis.cc
  bioMultiscaleRVE.cc
  bioRVEUtil.cc
  bioUtil.cc
  lasCSR.cc
  lasCSRBuilder.cc
  lasSparskit.cc
  lasSparskitExterns.cc)

find_package(sparskit REQUIRED)
find_package(scorecutil REQUIRED)
find_package(amsi COMPONENTS util analysis multiscale REQUIRED)

add_library(micro_fo ${micro_fo_src_files} ${micro_fo_public_headers})
set_property(TARGET micro_fo PROPERTY VERSION ${BIOTISSUE_VERSION})

target_include_directories(micro_fo PUBLIC
                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                    $<INSTALL_INTERFACE:include/biotissue>
                                    ${SCORECUTIL_INCLUDE_DIR})

target_link_libraries(micro_fo PUBLIC
                      ${SPARSKIT_LIBRARY}
                      ${SCORECUTIL_LIBRARIES}
                      amsi::util amsi::analysis amsi::multiscale)

install(TARGETS micro_fo
        EXPORT biotissueTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include/biotissue)

install(FILES ${micro_fo_header_files} DESTINATION include/biotissue)

