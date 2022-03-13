cmake_minimum_required(VERSION 3.21.1)
project(MuMFiM)
include(CTest)
set(BIOTISSUE_MAJOR_VERSION 0)
set(BIOTISSUE_MINOR_VERSION 2)
set(BIOTISSUE_PATCH_VERSION 7)
set(BIOTISSUE_VERSION ${BIOTISSUE_MAJOR_VERSION}.${BIOTISSUE_MINOR_VERSION}.${BIOTISSUE_PATCH_VERSION})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_DEBUG_POSTFIX d)
if(LOGRUN)
  add_definitions(-DLOGRUN)
endif()

# Force the project to use the CXX14 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add warning and error flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")

if(BUILD_TESTS)
  message(STATUS "BUILDING TESTS")
  function(add_mpi_test name file no_mpi_proc)
    message(STATUS "Adding test " ${name})
    add_executable(${name} ${file})
    set(test_parameters ${CMAKE_CURRENT_BINARY_DIR}/${name} ${ARGN})
    add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${no_mpi_proc}
  ${MPIEXEC_PREFLAGS} ${test_parameters} )
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

option(BUILD_EXTERNAL "build external libs with fetch content (except scorec)" ON)
if(${BUILD_EXTERNAL})
  add_subdirectory(external)
endif()

add_subdirectory(src)
if(${BUILD_TESTING})
  add_subdirectory(test)
endif()
