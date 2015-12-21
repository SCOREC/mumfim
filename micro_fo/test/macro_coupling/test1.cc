#include "MacroCoupling.h"
#include "RVE.h"
#include <mpi.h>
#include <PCU.h>
#include <iostream>
#include <vector>
int main(int argc, char * argv[])
{
  int result = 0;
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  std::vector<apf::Vector3> tet;
  tet.push_back(apf::Vector3(-1.22474487139159 ,  0.0,  0.0));
  tet.push_back(apf::Vector3( 0.408248290463863,  0.0, -1.15470053837925));
  tet.push_back(apf::Vector3( 0.408248290463863,  1.0,  0.577359269189626));
  tet.push_back(apf::Vector3( 0.408248290463863, -1.0,  0.577359269189626));
  bio::MacroInfo cplng(bio::makeSingleEntityMesh(apf::Mesh::TET,&tet[0]));
  apf::DynamicVector du_fe(12);
  du_fe.zero();
  du_fe(3) = du_fe(6) = du_fe(9) = 0.1;
  bio::RVE rve(3);
  apf::DynamicMatrix drve_dfe;
  cplng.calcdRVEdFE(drve_dfe,&rve);
  apf::DynamicVector du_rve;
  apf::multiply(drve_dfe,du_fe,du_rve);
  PCU_Comm_Free();
  MPI_Finalize();
  return result;
}
