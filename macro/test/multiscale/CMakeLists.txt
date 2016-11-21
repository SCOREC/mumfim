add_executable(multiscale main.cc)
find_package(amsi COMPONENTS util meta interface control)
target_link_libraries(multiscale
                      biotissue
                      micro_fm
                      micro_fo
                      util
                      meta
                      interface
                      control)
set_target_properties(multiscale PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
install(TARGETS multiscale
        RUNTIME DESTINATION bin)
