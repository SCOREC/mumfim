#include "Analysis.h"
#include "FiberNetwork.h"
#include "RVE.h"
#include <mpi.h>
#include <cassert>
#include <fstream>
#include <iostream>
int main(int argc, char * argv[])
{
  assert(argv[1]);
  int result = 0;
  MPI_Init(&argc,&argv);
  // bio::FiberNetwork * fn = bio::loadFromFile(std::string(argv[1]));
  // bio::RVE rve;
  // bio::FiberRVEAnalysis alsys(fn,&rve);
  // alsys.init();
  // alsys.run(1);
  MPI_Finalize();
  return result;
}

    
