#include "FiberNetworkIO.h"
#include <apf.h>
#include <PCU.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <mpi.h>
int main(int argc, char * argv[])
{
  assert(argv[1]);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
//! [load from file stream]
  std::ifstream strm(argv[1]);
  apf::Mesh2 * network = bio::NetworkLoader().fromStream(strm);
//! [load from file stream]
  strm.close();
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
