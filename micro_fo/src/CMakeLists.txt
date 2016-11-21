set(micro_fo_public_headers
  RepresentVolElem.h
  FiberNetwork.h
  RVE_Util.h
  globals.h
  MicroFOMultiscaleTypes.h
  SparseMatrix.h)

set(micro_fo_src_files
  RepresentVolElem.cc
  FiberNetwork.cc
  FiberNetworkIO.cc
  AlignFiberNetwork.cc
  main_solver_psi.cc
  solver_gmres_df.cc
  solver_beam.cc
  create_element_new.cc
  make_structure_fict_i2.cc
  avgvolstress.cc
  RVE_Util.cc
  SparseMatrix.cc
  Sparskit_Externs.cc)

find_package(sparskit REQUIRED)
find_package(scorecutil REQUIRED)
find_package(amsi COMPONENTS util meta interface control)

add_library(micro_fo ${micro_fo_src_files} ${micro_fo_public_headers})
set_property(TARGET micro_fo PROPERTY VERSION ${BIOTISSUE_VERSION})

target_include_directories(micro_fo PUBLIC
                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                    $<INSTALL_INTERFACE:include/biotissue>
                                    ${SCORECUTIL_INCLUDE_DIR})

target_link_libraries(micro_fo PUBLIC
                      ${SPARSKIT_LIBRARY}
                      ${SCORECUTIL_LIBRARIES}
                      util meta control interface)

install(TARGETS micro_fo
        EXPORT biotissueTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include/biotissue)

install(FILES ${micro_fo_header_files} DESTINATION include/biotissue)

