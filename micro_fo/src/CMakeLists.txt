set(micro_fo_public_headers
  FiberNetwork.h
  FiberReactions.h
  RepresentVolElem.h
  RVE_Util.h
  globals.h
  MicroFOMultiscaleTypes.h 
  SparseMatrix.h)

set(micro_fo_src_files
  avgvolstress.cc
  FiberNetwork.cc
  FiberNetworkIO.cc
  FiberReactions.cc
  RepresentVolElem.cc
  main_solver_psi.cc
  solver_gmres_df.cc
  solver_beam.cc
  create_element_new.cc
  make_structure_fict_i2.cc
  RVE_Util.cc 
  SparseMatrix.cc
  Sparskit_Externs.cc)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

if(NOT (TARGET sparskit))
  find_package(sparskit REQUIRED)
endif()
find_package(amsi)

add_library(micro_fo ${micro_fo_src_files} ${micro_fo_public_headers})
set_property(TARGET micro_fo PROPERTY VERSION ${BIOTISSUE_VERSION})

target_include_directories(micro_fo PUBLIC
				    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
				    $<INSTALL_INTERFACE:include/biotissue>)

target_link_libraries(micro_fo PUBLIC
		      ${SPARSKIT_LIBRARIES}
		      ${SCORECUTIL_LIBRARIES}
		      amsi_util 
		      amsi_meta 
		      amsi_control 
		      amsi_interface)

install(TARGETS micro_fo
	EXPORT biotissueTargets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include/biotissue)

install(FILES ${micro_fo_header_files} DESTINATION include/biotissue)

