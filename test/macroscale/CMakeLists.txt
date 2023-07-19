if (NOT ${BUILD_EXTERNAL})
    find_package(Catch2 REQUIRED)
endif ()

include(Catch)
add_executable(macro_unit_tests test_macro_main.cc NeohookeanMaterial.cc MultiscaleMaterial.cc UpdatedLagrangianMaterialIntegrator.cc ApfMatrixUtils.cc TestDeformations.cc TestSupport.cc TestSupport.h TorchMaterial.cc)
target_link_libraries(macro_unit_tests PUBLIC Catch2::Catch2 biotissue)
catch_discover_tests(macro_unit_tests)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
