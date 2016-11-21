add_executable(singlescale main.cc)
find_package(amsi COMPONENTS util meta interface control)
target_link_libraries(singlescale
                      biotissue
                      util
                      meta
                      interface
                      control)
set_target_properties(singlescale PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
install(TARGETS singlescale
        RUNTIME DESTINATION bin)
