find_package(amsi COMPONENTS util analysis multiscale)
add_executable(bio_multiscale main.cc)
target_link_libraries(bio_multiscale biotissue micro_fm micro_fo util analysis multiscale)
set_target_properties(bio_multiscale PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
install(TARGETS bio_multiscale
        RUNTIME DESTINATION bin)
