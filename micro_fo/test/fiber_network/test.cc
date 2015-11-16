#include "FiberNetworkIO.h"
#include "FiberNetwork.h"

#include <apf.h>
#include <PCU.h>

#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <mpi.h>

int main(int argc, char * argv[])
{
  assert(argv[1]);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();

  std::ifstream strm(argv[1]);
  apf::Mesh * ntwk = bio::NetworkLoader().fromStream(strm);
  strm.close();

  bio::FiberNetwork fn(ntwk);

  std::cout << "Network is " << fn.getDim() << "d" << std::endl;
  
  std::vector<double> lngths;
  fn.calcFiberLengths(lngths);
  std::cout << "Lengths for " << lngths.size() << " fibers" << std::endl;

  double ttl_lngth = std::accumulate(lngths.begin(),lngths.end(),0.0);
  std::cout << "Total length: " << ttl_lngth << std::endl;
  std::cout << "Avg length: " << ttl_lngth / (double)lngths.size() << std::endl;

  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
