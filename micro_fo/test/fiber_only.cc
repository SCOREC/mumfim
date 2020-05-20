#include <amsiAnalysis.h>
#include <amsiDetectOscillation.h>
#include <apf.h>
#include <lasCSRCore.h>
#include <lasCorePETSc.h>
#include <mpi.h>
#include <iostream>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOParams.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "bioVerbosity.h"
#include "bioFiberNetworkLibrary.h"
#include <Kokkos_Core.hpp>
void stressToMat(double * stress_arr, apf::Matrix3x3 & stress)
{
  stress[0][0] = stress_arr[0];
  stress[1][1] = stress_arr[1];
  stress[2][2] = stress_arr[2];
  stress[1][2] = stress_arr[3];
  stress[2][1] = stress_arr[3];
  stress[0][2] = stress_arr[4];
  stress[2][0] = stress_arr[4];
  stress[0][1] = stress_arr[5];
  stress[1][0] = stress_arr[5];
}
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
#ifdef MICRO_USING_PETSC
  las::initPETScLAS(&argc, &argv, MPI_COMM_WORLD);
#endif
  Kokkos::ScopeGuard kokkos(argc, argv);
  if(argc != 2)
  {
    std::cerr<<"Usage: "<<argv[0]<<" job.yaml"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  std::vector<bio::MicroCase> cases;
  bio::loadMicroFOFromYamlFile(argv[1], cases);
  bio::printMicroFOCase(cases[0]);
  std::string file_name = cases[0].pd.meshFile;
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bio::FiberNetworkLibrary network_library;
  network_library.load(file_name,file_name+".params",0,0);
  auto fiber_network = network_library.getOriginalNetwork(0,0);
  // I'm not confident that the move thing here works as intended
  auto an = bio::createFiberRVEAnalysis(std::move(fiber_network), std::move(cases[0].ss));

  double stress[6];
  double C[36];
  apf::Matrix3x3 strss;
  Kokkos::Timer timer;
  double time = timer.seconds();
  timer.reset();
  bool result = an->run(cases[0].pd.deformationGradient, stress);
  an->computeMaterialStiffness(C);
  std::cout<<"Stress from run"<<std::endl;
  stressToMat(stress, strss);
  std::cout<<strss<<std::endl;
  for(int i=0; i<6; ++i)
  {
    for(int j=0; j<6; ++j)
    {
      std::cout<<C[i*6+j]<<" ";
    }
    std::cout<<"\n";
  }
  std::cout<<std::endl;
  // DEBUG
  bio::DeformationGradient dg(1,0,0,0,1,0,0,0,1);
  result = an->run(dg, stress);
  an->computeMaterialStiffness(C);
  std::cout<<"Stress from run"<<std::endl;
  stressToMat(stress, strss);
  std::cout<<strss<<std::endl;
  for(int i=0; i<6; ++i)
  {
    for(int j=0; j<6; ++j)
    {
      std::cout<<C[i*6+j]<<" ";
    }
    std::cout<<"\n";
  }
  std::cout<<std::endl;
  // end DEBUG
  double time2 = timer.seconds();
  std::cout<<"Took: "<<time2 << " seconds."<<std::endl;
  if (!result)
  {
    std::cerr << "The microscale analysis failed to converge" << std::endl;
  }
  std::stringstream sout;
  sout << "rnk_" << rank << "_fn_" << an->getFn()->getRVEType();
  apf::writeVtkFiles(sout.str().c_str(), an->getFn()->getNetworkMesh(), 1);
#ifdef MICRO_USING_PETSC
  las::finalizePETScLAS();
#endif
  amsi::freeAnalysis();
  return 0;
}
