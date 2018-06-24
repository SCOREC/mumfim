set(biotissue_public_headers
  bioAnalysis.h
  bioAnalysisIO.h
  bioAnalysisIO_impl.h
  bioBatheMultiscale.h
  bioHolmesMowIntegrator.h
  bioLinearTissue.h
  bioMultiscaleConvergence.h
  bioMultiscaleIntegrator.h
  bioMultiscaleTissue.h
  bioNeoHookeanIntegrator.h
  bioNonlinearTissue.h
  bioStiffnessVariation.h
  bioStiffnessVariation_impl.h
  bioTissueAnalysis.h
  bioTissueAnalysis.h
  bioTissueAnalysis_impl.h
  bioTissueMultiscaleAnalysis.h
  bioULMultiscaleHydrostaticPressureIntegrator.h
  bioULMultiscaleIntegrator.h
  bioVariableRecovery.h
  bioVariableRecovery_impl.h
  bioVolumeConstraint.h
  bioVolumeConvergence.h)

set(biotissue_src_files
  bioAnalysis.cc
  bioLinearTissue.cc
  bioMultiscaleTissue.cc
  bioNonlinearTissue.cc
  bioStiffnessVariation.cc
  bioTissueAnalysis.cc
  bioTissueMultiscaleAnalysis.cc
  bioVolumeConstraint.cc
  bioVolumeConvergence.cc)

find_package(amsi COMPONENTS util analysis multiscale)

add_library(biotissue ${biotissue_src_files} ${biotissue_public_headers})
target_link_libraries(biotissue
                      micro_fo
                      amsi::util
                      amsi::analysis
                      amsi::multiscale)

target_include_directories(biotissue PUBLIC
                                     "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                                     "$<INSTALL_INTERFACE:include/biotissue>")

install(TARGETS biotissue DESTINATION lib)
install(FILES ${biotissue_public_headers} DESTINATION include/biotissue)

