#include "Analysis.h"
#include "FiberNetwork.h"
#include "RVE.h"
#include <mpi.h>
#include <PCU.h>
#include <cassert>
#include <fstream>
#include <iostream>
int main(int argc, char * argv[])
{
  assert(argc == 2);
  int result = 0;
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  std::string fn_nm(argv[1]);
  bio::FiberRVEAnalysis * rve = bio::makeAnalysis(fn_nm);
  apf::DynamicVector du(24);
  du.zero();
  du[0] = du[6] = du[12] = du[18] = -0.0306186217847201;
  du[3] = du[9] = du[15] = du[21] =  0.0306186217847201;
  bio::displaceRVE(rve->rve,du);
  bio::forwardRVEDisplacement(rve->rve,rve->fn);
  bio::FiberRVEIteration itr(rve);
  bio::FiberRVEConvergence cnv(rve,1e-8);
  num::numericalSolve(&itr,&cnv);
  bio::destroyAnalysis(rve);
  PCU_Comm_Free();
  MPI_Finalize();
  return result;
}

    
