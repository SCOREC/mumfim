find_package(amsi COMPONENTS util analysis multiscale REQUIRED)
find_package(las)
#find_package(Kokkos CONFIG REQUIRED)
add_executable(bio_multiscale main.cc)
target_link_libraries(bio_multiscale biotissue micro_fo las amsi::util amsi::analysis amsi::multiscale)
#target_include_directories(bio_multiscale ${Kokkos_INCLUDE_DIRS})
set_target_properties(bio_multiscale PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
install(TARGETS bio_multiscale
        RUNTIME DESTINATION bin)
