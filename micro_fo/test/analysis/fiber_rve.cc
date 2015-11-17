#include "Solver.h"
#include "FiberNetworkIO.h"
#include "FiberNetwork.h"
#include "RVE.h"

#include <cassert>
#include <fstream>
#include <iostream>

int main(int argc, char * argv[])
{
  assert(argv[1]);
  std::ifstream strm(argv[1]);
  apf::Mesh2 * network = bio::NetworkLoader().fromStream(strm);
  bio::FiberNetwork fn(network);
  bio::RVE rve;
  bio::FiberRVEAnalysis alsys(&fn,&rve);
  alsys.run(1);
}

    
