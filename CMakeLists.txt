cmake_minimum_required(VERSION 3.21.1)
project(MuMFiM)
set(BIOTISSUE_MAJOR_VERSION 0)
set(BIOTISSUE_MINOR_VERSION 2)
set(BIOTISSUE_PATCH_VERSION 7)
set(BIOTISSUE_VERSION ${BIOTISSUE_MAJOR_VERSION}.${BIOTISSUE_MINOR_VERSION}.${BIOTISSUE_PATCH_VERSION})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_DEBUG_POSTFIX d)
if(LOGRUN)
  add_definitions(-DLOGRUN)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)



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

option(MUMFIM_ENABLE_ADIOS2 "enable adios2 test cases" OFF)

add_subdirectory(src)
if(BUILD_TESTS)
  include(CTest)
  enable_testing()
  include(add-mpi-test)
  add_subdirectory(test)
endif(BUILD_TESTS)

