# Once done this will define
#  SPARSKIT_FOUND - System has SPARSKIT
#  SPARSKIT_LIBRARIES - The libraries needed to use SPARSKIT
#  SPARSKIT_DEFINITIONS - Compiler switches required for using SPARSKIT

checkSetParam(SPARSKIT_DIR FALSE)

find_library(SPARSKIT_LIBRARY Sparskit
             HINTS ${SPARSKIT_DIR}
             PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPARSKIT DEFAULT_MSG SPARSKIT_LIBRARY)
mark_as_advanced( SPARSKIT_LIBRARY )
