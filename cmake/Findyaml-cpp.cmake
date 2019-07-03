#include(util)

find_package(PkgConfig REQUIRED QUIET)

if(yaml-cpp_FIND_REQUIRED)
  set(_yaml-cpp_OPTS "REQUIRED")
endif()
if(yaml-cpp_FIND_QUIETLY)
  set(_yaml-cpp_OPTS "QUIET")
endif()
if(yayaml-cpp_FIND_REQUIRED AND yaml-cpp_FIND_QUIETLY)
  set(_yaml-cpp_OPTS "REQUIRED QUIET")
endif()

set(ENV{PKG_CONFIG_PATH}  "${yaml-cpp_DIR}")

pkg_check_modules(YAML-CPP ${_yaml-cpp_OPTS} yaml-cpp)
#message(INFO ${yaml-cpp_INCLUDEDIR})
#message(INFO ${yaml-cpp_INCLUDE_DIRS})

if(YAML-CPP_FOUND)
  set(yaml-cpp_INCLUDE_DIRS "${YAML-CPP_INCLUDEDIR}")
  #set(yaml-cpp_INCLUDE_DIRS "${YAML-CPP_INCLUDE_DIRS}")
  #set(yaml-cpp_LIBRARIES "${YAML-CPP_STATIC_LDFLAGS};${YAML-CPP_STATIC_LIBRARIES}")
  set(yaml-cpp_LIBRARIES "${YAML-CPP_LDFLAGS};${YAML-CPP_LIBRARIES}")
endif(YAML-CPP_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(yaml-cpp DEFAULT_MSG yaml-cpp_LIBRARIES yaml-cpp_INCLUDE_DIRS)
mark_as_advanced(yaml-cpp_LIBRARIES yaml-cpp_INCLUDE_DIRS)
