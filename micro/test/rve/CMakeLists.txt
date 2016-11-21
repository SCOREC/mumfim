find_package(amsi REQUIRED)

add_executable(rve rve.cc)
target_link_libraries(rve micro_fo las num fea ${CORE_LIBRARIES})