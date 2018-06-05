#add_mpi_test(convert_stress_derivs "convert_stress_derivs.cc;io.cc" 1)
#target_link_libraries(convert_stress_derivs micro_fo)
#configure_file(ds_dx_rve.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(strs.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(dV_dx_rve.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(vol_cnv.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(dS_dx_rve.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

add_mpi_test(drvedfe "dRVEdFE.cc;io.cc" 1)
target_link_libraries(drvedfe micro_fo)
configure_file(dx_rve_dx_fe1.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(dx_rve_dx_fe2.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(fn.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# depends on fn.test
add_mpi_test(dRdxrve "dR_dx_rve.cc;io.cc" 1)
target_link_libraries(dRdxrve micro_fo)
configure_file(dR_dx_rve.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

add_mpi_test(dx_fn_dx_rve "dx_fn_dx_rve.cc;io.cc" 1)
target_link_libraries(dx_fn_dx_rve micro_fo)
configure_file(K_dR_dx_rve.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(dx_fn_dx_rve.test ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
