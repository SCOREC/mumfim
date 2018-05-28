add_mpi_test(convert_stress_derivs "convert_stress_derivs.cc;io.cc" 1)
target_link_libraries(convert_stress_derivs micro_fo)

add_mpi_test(drvedfe "dRVEdFE.cc;io.cc" 1)
target_link_libraries(drvedfe micro_fo)
