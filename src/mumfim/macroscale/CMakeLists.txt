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
        VolumeConvergence.cc
        UpdatedLagrangianMaterialIntegrator.cc
        UpdatedLagrangianMaterialIntegrator.h
        materials/NeoHookean.cc
        ApfMatrixUtils.cc
        materials/MultiscaleMaterial.cc
        PetscSNES.cpp)

if (MUMFIM_ENABLE_Torch)
    find_package(Torch REQUIRED)
    list(APPEND biotissue_src_files materials/TorchMaterial.cc)
endif ()
if (NOT ${BUILD_EXTERNAL})
    find_package(amsi COMPONENTS util analysis multiscale)
endif ()

add_library(biotissue ${biotissue_src_files})

target_link_libraries(biotissue
        PUBLIC
        micro_fo
        amsi::util
        amsi::analysis
        amsi::multiscale
        amsi::amsi)

if (MUMFIM_ENABLE_Torch)
    target_link_libraries(biotissue PUBLIC ${TORCH_LIBRARIES})
    target_compile_definitions(biotissue PUBLIC -DMUMFIM_ENABLE_TORCH)
endif ()

target_include_directories(biotissue PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        "$<INSTALL_INTERFACE:include>")

target_compile_features(biotissue PUBLIC cxx_std_17)
install(TARGETS biotissue DESTINATION lib)
