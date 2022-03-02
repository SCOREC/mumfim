if(NOT ${BUILD_EXTERNAL})
find_package(amsi COMPONENTS multiscale REQUIRED)
endif()

# This test should be moved to las! not sure why bio specific stuff is used here
# we don't actually check anything in this test anyways!
#add_mpi_test(csr test_csr.cc 1)
#target_link_libraries(csr micro_fo)

add_mpi_test(rve test_rve.cc 1)
target_link_libraries(rve micro_fo)

add_mpi_test(fiber_only fiber_only.cc 1 ${CMAKE_CURRENT_BINARY_DIR}/fiber_only.yaml)
target_link_libraries(fiber_only micro_fo)

add_mpi_test(batched_fiber_only batched_fiber_only.cc 1 ${CMAKE_CURRENT_BINARY_DIR}/fiber_only.yaml)
target_link_libraries(batched_fiber_only PUBLIC micro_fo) 

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/fiber_only.yaml" DESTINATION
          "${CMAKE_CURRENT_BINARY_DIR}")

add_mpi_test(test_global_stiffness ./test_global_stiffness.cc 1)
target_link_libraries(test_global_stiffness micro_fo)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/test_global_stiffness_data" DESTINATION
          "${CMAKE_CURRENT_BINARY_DIR}/")

# FIXME This test doesn't actually check anything...need to update this test case to make sure
# that we actually test that the elemental stiffness matrix is getting computed properly!
add_mpi_test(test_elemental_stiffness ./test_elemental_stiffness.cc 1 ${CMAKE_CURRENT_BINARY_DIR}/fiber_only.yaml)
target_link_libraries(test_elemental_stiffness micro_fo)

#add_mpi_test(test_mass_integrator ./test_mass_integrator.cc 1)
#target_link_libraries(test_mass_integrator micro_fo)

file(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/fiber_networks
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


# new unit tests require catch 2
if(NOT ${BUILD_EXTERNAL})
find_package(Catch2)
endif()
if(catch2_FOUND)
include(ParseAndAddCatchTests)

set(TEST_SOURCES test_main.cc
                 test_reorder_mesh.cc
                 test_kokkos_view_of_views.cc
   )
add_executable(unit_tests ${TEST_SOURCES})
target_link_libraries(unit_tests PUBLIC Catch2::Catch2 micro_fo)
ParseAndAddCatchTests(unit_tests)
endif()
