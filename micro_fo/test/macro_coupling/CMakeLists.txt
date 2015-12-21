find_package(core COMPONENTS sim)

add_executable(macro_coupling test1.cc)
target_link_libraries(macro_coupling micro_fo ${CORE_LIBRARIES})