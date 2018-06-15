#include "bioRVE.h"
#include <PCU.h>
#include <mpi.h>
int main(int argc, char * argv[])
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  bio::RVE rve(3);
  apf::DynamicVector xyz_rve_0;
  bio::getRVEReferenceCoords(&rve,xyz_rve_0);
  apf::DynamicVector u_rve;
  bio::getRVEDisplacement(&rve,u_rve);
  PCU_Comm_Free();
  MPI_Finalize();
}
