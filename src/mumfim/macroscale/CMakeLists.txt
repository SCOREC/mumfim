set(biotissue_src_files
        Analysis.cc
        LinearTissue.cc
        MultiscaleTissue.cc
        ModelTraits.cc
        NonlinearTissue.cc
        ReadStochasticField.cc
        RVECoupling.cc
        StiffnessVariation.cc
        TissueAnalysis.cc
        MultiscaleTissueAnalysis.cc
        SinglescaleTissueAnalysis.cc
        VolumeConstraint.cc
        VolumeConvergence.cc)

if (NOT ${BUILD_EXTERNAL})
    find_package(amsi COMPONENTS util analysis multiscale)
endif ()

add_library(biotissue ${biotissue_src_files})
#target_link_libraries(biotissue
#        micro_fo
#        amsi::util
#        amsi::analysis
#        amsi::multiscale)
target_link_libraries(biotissue
        PUBLIC
        micro_fo
        amsi::util
        amsi::analysis
        amsi::multiscale
        amsi::amsi)

target_include_directories(biotissue PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        "$<INSTALL_INTERFACE:include>")

install(TARGETS biotissue DESTINATION lib)
