#include "RVE.h"
#include <mpi.h>
#include <PCU.h>
int main(int argc, char * argv[])
{
  assert(argv[1]);
  int result = 0;
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  bio::RVE rve;
  int ndes = rve.numNodes();
  PCU_Comm_Free();
  MPI_Finalize();
  return result;
}
