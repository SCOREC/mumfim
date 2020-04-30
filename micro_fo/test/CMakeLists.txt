find_package(amsi COMPONENTS multiscale REQUIRED)

add_mpi_test(csr test_csr.cc 1)
target_link_libraries(csr micro_fo)

add_mpi_test(rve test_rve.cc 1)
target_link_libraries(rve micro_fo)

add_mpi_test(fiber_only fiber_only.cc 1)
target_link_libraries(fiber_only micro_fo)

add_mpi_test(test_global_stiffness ./test_global_stiffness.cc 1)
target_link_libraries(test_global_stiffness micro_fo)
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/test_global_stiffness_data" DESTINATION
          "${CMAKE_CURRENT_BINARY_DIR}")
#add_custom_command(TARGET test_global_stiffness POST_BUILD
#                   COMMAND ${CMAKE_COMMAND} -E copy_directory
#                   "${CMAKE_CURRENT_SOURCE_DIR}/test_global_stiffness"
#                   "${CMAKE_CURRENT_BINARY_DIR}/test_global_stiffness")
#                   #"$<TARGET_FILE_DIR:test_global_stiffness>/test_global_stiffness")

add_mpi_test(test_elemental_stiffness ./test_elemental_stiffness.cc 1)
target_link_libraries(test_elemental_stiffness micro_fo)

#add_mpi_test(test_mass_integrator ./test_mass_integrator.cc 1)
#target_link_libraries(test_mass_integrator micro_fo)

file(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/fiber_networks
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


# new unit tests require catch 2
find_package(Catch2 REQUIRED)
include(ParseAndAddCatchTests)

set(TEST_SOURCES test_main.cc
                 test_kokkos_view_of_views.cc
   )
add_executable(unit_tests ${TEST_SOURCES})
target_link_libraries(unit_tests PUBLIC Catch2::Catch2 micro_fo)
ParseAndAddCatchTests(unit_tests)
