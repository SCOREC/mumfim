add_subdirectory(mumfim/macroscale)
add_subdirectory(mumfim/microscale)

if (NOT ${BUILD_EXTERNAL})
    find_package(amsi COMPONENTS util analysis multiscale REQUIRED)
    find_package(ModelTraits)
    find_package(las)
endif ()

add_executable(singlescale_analysis singlescale_analysis.cc)
target_link_libraries(singlescale_analysis biotissue amsi::util amsi::analysis ModelTraits::ModelTraits ModelTraits::IO)
target_compile_features(singlescale_analysis PUBLIC cxx_std_20)
install(TARGETS singlescale_analysis RUNTIME DESTINATION bin)

add_executable(multiscale_analysis multiscale_analysis.cc)
target_compile_features(multiscale_analysis PUBLIC cxx_std_20)
target_link_libraries(multiscale_analysis biotissue micro_fo amsi::util amsi::analysis amsi::multiscale ModelTraits::ModelTraits ModelTraits::IO)
install(TARGETS multiscale_analysis RUNTIME DESTINATION bin)
