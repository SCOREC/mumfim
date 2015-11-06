set(micro_fo_public_headers
  FiberNetworkIO.h
  FiberReactions.h
  MacroCoupling.h
  SparseStructure.h)

set(micro_fo_src_files
  dRVEdFE.cc
  FiberNetworkIO.cc
  FiberReactions.cc
  SparseStructure.cc
  Sparskit_Externs.cc)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

if(NOT (TARGET sparskit))
  find_package(sparskit REQUIRED)
endif()
#find_package(amsi)
find_package(CoreSim REQUIRED)

add_library(micro_fo ${micro_fo_src_files} ${micro_fo_public_headers})
set_property(TARGET micro_fo PROPERTY VERSION ${BIOTISSUE_VERSION})

target_include_directories(micro_fo PUBLIC
				    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
				    $<INSTALL_INTERFACE:include/biotissue>
                                    ${CORE_INCLUDE_DIRS})

target_link_libraries(micro_fo PUBLIC
		      ${SPARSKIT_LIBRARIES}
		      ${SCORECUTIL_LIBRARIES}
		      ${CORE_LIBRARIES}
                      ${MPI_LIBRARIES})

install(TARGETS micro_fo
	EXPORT biotissueTargets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include/biotissue)

install(FILES ${micro_fo_header_files} DESTINATION include/biotissue)

