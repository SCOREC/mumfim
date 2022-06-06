#define CATCH_CONFIG_RUNNER
#include <mpi.h>
#include "catch2/catch.hpp"
#include <Kokkos_Core.hpp>
#include <PCU.h>
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  auto result= Catch::Session().run(argc, argv);
  PCU_Comm_Free();
  MPI_Finalize();
  return result;
}
