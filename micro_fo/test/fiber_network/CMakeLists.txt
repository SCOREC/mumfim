#add_executable(converter converter.cc)
#target_link_library(converter micro_fo)


add_executable(load load.cc)
target_link_libraries(load micro_fo ${CORE_LIBRARIES})

add_executable(fiber_network test.cc)
target_link_libraries(fiber_network micro_fo ${CORE_LIBRARIES})