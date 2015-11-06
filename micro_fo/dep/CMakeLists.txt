enable_language(Fortran)
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")

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
add_library(sparskit ${SOURCE_FILES})

# For install, put the created AMSI library in installpath/lib
install (TARGETS sparskit DESTINATION lib)
