#include "mumfim/microscale/RVE.h"
#include <PCU.h>
#include <mpi.h>
int main(int argc, char * argv[])
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  mumfim::RVE rve(3);
  apf::DynamicVector xyz_rve_0;
  mumfim::getRVEReferenceCoords(&rve,xyz_rve_0);
  apf::DynamicVector u_rve;
  mumfim::getRVEDisplacement(&rve,u_rve);
  PCU_Comm_Free();
  MPI_Finalize();
}
