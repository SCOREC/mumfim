cmake_minimum_required(VERSION 3.0)
project(biotissue)
set(BIOTISSUE_MAJOR_VERSION 1)
set(BIOTISSUE_MINOR_VERSION 0)
set(BIOTISSUE_PATCH_VERSION 0)
set(BIOTISSUE_VERSION ${BIOTISSUE_MAJOR_VERSION}.${BIOTISSUE_MINOR_VERSION}.${BIOTISSUE_PATCH_VERSION})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_FIND_LIBRARY_SUFFIX ".a")

if(LOGRUN)
  add_definitions(-DLOGRUN)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  list(APPEND CMAKE_CXX_FLAGS "-qlanglvl=extended")
endif()

add_subdirectory(micro_fo)
#add_subdirectory(micro_fm)
add_subdirectory(macro)
