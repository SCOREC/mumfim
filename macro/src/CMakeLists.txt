set(biotissue_public_headers
  bioAnalysis.h
  bioAnalysisIO.h
  bioAnalysisIO_impl.h
  bioMultiscaleTissue.h
  bioTissueAnalysis.h
  bioMultiscaleIntegrator.h
  bioULMultiscaleIntegrator.h
  bioNeoHookeanIntegrator.h
  bioHolmesMowIntegrator.h
  bioBatheMultiscale.h
  bioPreserveVolConstraintIntegrator.h
  bioVolumeConstraintIntegrator.h
  bioLinearTissue.h
  bioNonlinearTissue.h
  bioTissueMultiscaleAnalysis.h
  bioVariableRecovery.h
  bioVariableRecovery_impl.h)

set(biotissue_src_files
  bioAnalysis.cc
  bioLinearTissue.cc
  bioMultiscaleTissue.cc
  bioNonlinearTissue.cc
  bioTissueMultiscaleAnalysis.cc)

find_package(amsi COMPONENTS util analysis multiscale)

add_library(biotissue ${biotissue_src_files} ${biotissue_public_headers})
target_link_libraries(biotissue
                      micro_fo
                      micro_fm
                      util
                      analysis
                      multiscale)

target_include_directories(biotissue PUBLIC
                                     "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                                     "$<INSTALL_INTERFACE:include/biotissue>")

install(TARGETS biotissue DESTINATION lib)
install(FILES ${biotissue_public_headers} DESTINATION include/biotissue)

