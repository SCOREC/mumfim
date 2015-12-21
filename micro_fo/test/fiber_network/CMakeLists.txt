#add_executable(converter converter.cc)
#target_link_library(converter micro_fo)

find_package(amsi REQUIRED)
find_package(Core COMPONENTS Sim)

add_executable(load load.cc)
target_link_libraries(load micro_fo ${CORE_LIBS})

add_executable(fiber_network test.cc)
target_link_libraries(fiber_network micro_fo las num fea ${CORE_LIBS})