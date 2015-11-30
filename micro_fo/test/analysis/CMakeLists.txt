find_package(amsi REQUIRED)

add_executable(fiber_rve fiber_rve.cc)
target_link_libraries(fiber_rve micro_fo las num fea ${CORE_LIBRARIES})