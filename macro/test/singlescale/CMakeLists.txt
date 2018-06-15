find_package(amsi COMPONENTS util analysis)
add_executable(bio_singlescale main.cc)
target_link_libraries(bio_singlescale biotissue amsi::util amsi::analysis)
set_target_properties(bio_singlescale PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
install(TARGETS bio_singlescale
        RUNTIME DESTINATION bin)
