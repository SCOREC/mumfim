#define CATCH_CONFIG_RUNNER
#include <mpi.h>
#include "catch2/catch.hpp"
#include <Kokkos_Core.hpp>
int main(int argc, char * argv[])
{
  int result;
  {
  Kokkos::ScopeGuard kokkos(argc, argv);
  result = Catch::Session().run(argc, argv);
  }
  return result;
}
