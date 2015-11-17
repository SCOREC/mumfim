#add_executable(vec vec.cc)
#target_link_libraries(vec micro_fo ${CORE_LIBRARIES})

add_executable(skls_init test.cc)
target_link_libraries(skls_init micro_fo ${CORE_LIBRARIES})