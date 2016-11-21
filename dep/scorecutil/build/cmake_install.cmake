# Install script for directory: /fasttmp/wtobin/develop/biotissue/dep/scorecutil

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/fasttmp/wtobin/develop/install/scorecutil/openmpi-1.10.0")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/build/libSCORECUtil.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/affine_space.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/BezierMappingBuilder.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/BezierMapping.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/GaussIntegrator.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/GaussLegendreSimplex.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/Integrator.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/GaussQuadrature.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/IntPt.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/LagrangeMappingBuilder.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/LagrangeMapping.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/MappingBuilder.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/Mapping.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/mPoint.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/MSMapping.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/mTensor2.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/mTensor2Symmetric.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/mTensor4.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/mVector.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/OctreeCreate2.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/Octree.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/o_internals.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/scorec_function_objects.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/scorecMD_memory.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/scorecSC_mem.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/scorecSListBase.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/scorecSSList.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/CountTime.h"
    "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/include/VolumeBuckets.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/module" TYPE FILE FILES "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/build/SCORECUTIL_1.0")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

file(WRITE "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/build/${CMAKE_INSTALL_MANIFEST}" "")
foreach(file ${CMAKE_INSTALL_MANIFEST_FILES})
  file(APPEND "/fasttmp/wtobin/develop/biotissue/dep/scorecutil/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
endforeach()
