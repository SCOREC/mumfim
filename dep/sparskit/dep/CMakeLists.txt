cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)
project (Sparskit)
set(PACKAGE_VERSION "1.0")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_Fortran_COMPILER mpif77)
else()
  set(CMAKE_Fortran_COMPILER mpif77)
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")

set(HEADER_FILES)

set(SOURCE_FILES 
  src/blas1.f
  src/ilu.f
  src/itaux.f
  src/iters.f
  src/ludcmp.f
  src/matvec.f
  src/pgmres.f
  src/pythag.f
  src/tql2.f
  src/tred2EISPACK.f
  src/utils.f )

# Add the library to create
add_library(Sparskit ${SOURCE_FILES} ${HEADER_FILES})

# For install, put the created AMSI library in installpath/lib
install (TARGETS Sparskit DESTINATION lib)

string(TOUPPER ${PROJECT_NAME} PACKAGE_NAME)
set(MODULEFILE ${PROJECT_NAME}_${PACKAGE_VERSION})
configure_file(moduletemplate ${MODULEFILE} @ONLY)

install(FILES build/${MODULEFILE} DESTINATION lib/module)

