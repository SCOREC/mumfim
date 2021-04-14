if (NOT ${BUILD_EXTERNAL})
    find_package(amsi COMPONENTS util analysis multiscale REQUIRED)
    find_package(las)
    find_package(ModelTraits)
endif ()
add_executable(bio_multiscale main.cc)

target_link_libraries(bio_multiscale biotissue micro_fo amsi::util amsi::analysis amsi::multiscale ModelTraits::ModelTraits ModelTraits::IO)
set_target_properties(bio_multiscale PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
install(TARGETS bio_multiscale
        RUNTIME DESTINATION bin)
