if (NOT ${BUILD_EXTERNAL})
    find_package(Catch2 REQUIRED)
endif ()

include(Catch)

set(TEST_SOURCES
    test_macro_main.cc
    NeohookeanMaterial.cc
    MultiscaleMaterial.cc
    UpdatedLagrangianMaterialIntegrator.cc
    ApfMatrixUtils.cc
    TestDeformations.cc
    TestSupport.cc
    TestSupport.h
  )

if(MUMFIM_ENABLE_Torch)
  list(APPEND TEST_SOURCES TorchMaterial.cc)
endif()


add_executable(macro_unit_tests ${TEST_SOURCES})
target_link_libraries(macro_unit_tests PUBLIC Catch2::Catch2 biotissue)
catch_discover_tests(macro_unit_tests)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
