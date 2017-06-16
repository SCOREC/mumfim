cmake_minimum_required(VERSION 3.0)
project(biotissue)
set(BIOTISSUE_MAJOR_VERSION 1)
set(BIOTISSUE_MINOR_VERSION 0)
set(BIOTISSUE_PATCH_VERSION 0)
set(BIOTISSUE_VERSION ${BIOTISSUE_MAJOR_VERSION}.${BIOTISSUE_MINOR_VERSION}.${BIOTISSUE_PATCH_VERSION})

list(APPEND CMAKE_MODULE_PATH $ENV{DEVROOT}/scripts/cmake ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(util)
set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so)
set(CMAKE_DEBUG_POSTFIX d)

if(LOGRUN)
  add_definitions(-DLOGRUN)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -Wall")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qflag=w -qlanglvl=extended0x")
  SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)
  SET_PROPERTY(GLOBAL PROPERTY TARGET_ARCHIVES_MAY_BE_SHARED_LIBS FALSE)
  set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
  set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
  set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif()

if(BUILD_TESTS)
  include(CTest)
  SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
  SET(GCC_COVERAGE_LINK_FLAGS "-lgcov")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  function(add_mpi_test name file no_mpi_proc)
    message(STATUS "Adding test " ${name})
    add_executable(${name} ${file})
    target_include_directories(${name} PUBLIC ${test_header_dir})
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      get_target_property(FLGS ${name} LINK_FLAGS_DEBUG)
      if(FLGS STREQUAL "FLGS-NOTFOUND")
        SET(FLGS "") # set to empty string
      endif()
      SET(FLGS "${FLGS} ${GCC_COVERAGE_LINK_FLAGS}" )
      set_target_properties(${name} PROPERTIES LINK_FLAGS_DEBUG ${FLGS} )
    endif()
    set(test_parameters -np ${no_mpi_proc} ${name} ${ARGN})
    add_test(NAME ${name} COMMAND "mpirun" ${test_parameters} )
  endfunction(add_mpi_test)
endif(BUILD_TESTS)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
                    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating API documentation with doxygen." VERBATIM)
endif()

find_package(SCOREC CONFIG REQUIRED)
add_subdirectory(micro_fo)
add_subdirectory(micro_fm)
add_subdirectory(macro)
