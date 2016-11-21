find_package(amsi REQUIRED)

add_executable(macro_coupling test1.cc)
target_link_libraries(macro_coupling micro_fo las num fea ${CORE_LIBRARIES})