if(NOT ${BUILD_EXTERNAL})
    find_package(Catch2 REQUIRED)
endif()

include(Catch)
add_executable(macro_unit_tests test_macro_main.cc NeohookeanMaterial.cc MultiscaleMaterial.cc UpdatedLagrangianMaterialIntegrator.cc ApfMatrixUtils.cc TestDeformations.cc TestSupport.cc TestSupport.h)
target_link_libraries(macro_unit_tests PUBLIC Catch2::Catch2 biotissue)
catch_discover_tests(macro_unit_tests)
