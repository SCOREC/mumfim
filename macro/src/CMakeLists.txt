set(biotissue_src_files
        bioAnalysis.cc
        bioLinearTissue.cc
        bioMultiscaleTissue.cc
        bioModelTraits.cc
        bioNonlinearTissue.cc
        bioReadStochasticField.cc
        bioRVECoupling.cc
        bioStiffnessVariation.cc
        bioTissueAnalysis.cc
        bioTissueMultiscaleAnalysis.cc
        bioVolumeConstraint.cc
        bioVolumeConvergence.cc)

if (NOT ${BUILD_EXTERNAL})
    find_package(amsi COMPONENTS util analysis multiscale)
endif ()

add_library(biotissue ${biotissue_src_files})
target_link_libraries(biotissue
        micro_fo
        amsi::util
        amsi::analysis
        amsi::multiscale)

target_include_directories(biotissue PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:include/biotissue>")

install(TARGETS biotissue DESTINATION lib)
