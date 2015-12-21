cmake_minimum_required(VERSION 3.0)
project(biotissue)

list(APPEND CMAKE_MODULE_PATH $ENV{DEVROOT}/scripts/cmake ${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_FIND_LIBRARY_SUFFIX .a .so)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  list(APPEND CMAKE_CXX_FLAGS "-qlanglvl=extended")
  set(CMAKE_EXE_LINKER_FLAGS -static)
endif()

add_subdirectory(micro_fo)
#add_subdirectory(micro_fm)
add_subdirectory(macro)
