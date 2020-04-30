#define CATCH_CONFIG_RUNNER
#include <mpi.h>
#include "catch2/catch.hpp"
//#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
//#endif
int main(int argc, char * argv[])
{
//#ifdef ENABLE_KOKKOS
  int result;
  Kokkos::initialize(argc, argv);
  {
//#endif
    result = Catch::Session().run(argc, argv);
//#ifdef ENABLE_KOKKOS
  }
  Kokkos::finalize();
//#endif
  return result;
}
