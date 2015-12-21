set(micro_fo_public_headers
  Analysis.h
  apfUtil.h
  dMdNi.h
  FiberNetworkIO.h
  FiberReactions.h
  FiberNetwork.h
  MacroCoupling.h
  RVE.h
  TrussIntegrator.h)

set(micro_fo_src_files
  Analysis.cc
  apfUtil.cc
  dRVEdFE.cc
  FiberNetworkIO.cc
  FiberNetwork.cc
  FiberReactions.cc
  MacroCoupling.cc
  RVE.cc)

find_package(amsi REQUIRED)
find_package(Core COMPONENTS Sim)

add_library(micro_fo ${micro_fo_src_files} ${micro_fo_public_headers})
set_property(TARGET micro_fo PROPERTY VERSION ${BIOTISSUE_VERSION})

target_include_directories(micro_fo PUBLIC
				    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
				    $<INSTALL_INTERFACE:include/biotissue>
                                    ${CORE_INCLUDE_DIRS})

target_link_libraries(micro_fo PUBLIC 
                               las
			       num
			       fea
                               ${CORE_LIBS})

install(TARGETS micro_fo
	EXPORT biotissueTargets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include/biotissue)

install(FILES ${micro_fo_header_files} DESTINATION include/biotissue)

