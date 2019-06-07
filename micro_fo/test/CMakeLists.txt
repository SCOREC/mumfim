find_package(amsi COMPONENTS multiscale REQUIRED)

add_mpi_test(csr csr.cc 1)
target_link_libraries(csr micro_fo)

add_mpi_test(rve testRVE.cc 1)
target_link_libraries(rve micro_fo)

add_mpi_test(fiber_only main.cc 1)
target_link_libraries(fiber_only micro_fo)

add_mpi_test(single_load_balance single_load_balance.cc 1)
target_link_libraries(single_load_balance micro_fo amsi::multiscale)

add_mpi_test(serialization serialization.cc 1)
target_link_libraries(serialization micro_fo)

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

#add_mpi_test(isolated_dof isolated_dof.cc 1)
#target_link_libraries(isolated_dof micro_fo)

#add_mpi_test(multi_disp multi_disp.cc 1)
#target_link_libraries(multi_disp micro_fo)

#add_mpi_test(singleDOF_disp singleDOF_disp.cc 1)
#target_link_libraries(singleDOF_disp micro_fo)

file(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/fiber_networks
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_mpi_test(interp_onto interp_onto.cc 1 )
target_link_libraries(interp_onto micro_fo)

add_subdirectory(coupling)
