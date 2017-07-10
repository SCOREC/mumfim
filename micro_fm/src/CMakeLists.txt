set(micro_fm_public_headers )

set(micro_fm_src_files )

find_package(amsi COMPONENTS util analysis multiscale REQUIRED)

add_library(micro_fm ${micro_fm_src_files} ${micro_fm_public_headers})

target_link_libraries(micro_fm amsi::util amsi::analysis amsi::multiscale)

target_include_directories(micro_fm PUBLIC
                                    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                                    "$<INSTALL_INTERFACE:include/biotissue>"
                                    ${CORE_INCLUDE_DIRS})

install(TARGETS micro_fm DESTINATION lib)

install(FILES ${micro_fm_public_headers} DESTINATION include/biotissue)

