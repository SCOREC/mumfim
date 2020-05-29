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
#include <Kokkos_DualView.hpp>
#include <bioBatchedFiberRVEAnalysisExplicit.h>
#include <memory>
template <typename T>
void stressToMat(int idx, T stress_view, apf::Matrix3x3 & stress)
{
  stress[0][0] = stress_view(idx,0);
  stress[1][1] = stress_view(idx,1);
  stress[2][2] = stress_view(idx,2);
  stress[1][2] = stress_view(idx,3);
  stress[2][1] = stress_view(idx,3);
  stress[0][2] = stress_view(idx,4);
  stress[2][0] = stress_view(idx,4);
  stress[0][1] = stress_view(idx,5);
  stress[1][0] = stress_view(idx,5);
}
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
#ifdef MICRO_USING_PETSC
  las::initPETScLAS(&argc, &argv, MPI_COMM_WORLD);
#endif
  Kokkos::ScopeGuard kokkos(argc, argv);
  if(argc != 3)
  {
    std::cerr<<"Usage: "<<argv[0]<<" job.yaml num_rves"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  auto BatchNum = std::atoi(argv[2]);
  std::vector<bio::MicroCase> cases;
  bio::loadMicroFOFromYamlFile(argv[1], cases);
  bio::printMicroFOCase(cases[0]);
  std::string file_name = cases[0].pd.meshFile;
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bio::FiberNetworkLibrary network_library;
  network_library.load(file_name,file_name+".params",0,0);
  std::vector<std::unique_ptr<bio::FiberNetwork>> fiber_networks;
  std::vector<std::shared_ptr<bio::MicroSolutionStrategy>> solution_strategies;
  std::shared_ptr<bio::MicroSolutionStrategy> shared_case{std::move(cases[0].ss)};

  for(int i=0; i<BatchNum; ++i)
  {
    fiber_networks.push_back(network_library.getCopy(0,0));
    solution_strategies.push_back(shared_case);
  }

  
  using ExeSpace = Kokkos::DefaultExecutionSpace;
  //using ExeSpace = Kokkos::Serial;
  //using Scalar = float;//bio::Scalar;
  using Scalar = bio::Scalar;
  bio::BatchedFiberRVEAnalysisExplicit<Scalar,bio::LocalOrdinal,ExeSpace> batched_analysis(std::move(fiber_networks),std::move(solution_strategies));

  Kokkos::DualView<Scalar*[3][3],ExeSpace> deformation_gradient("deformation gradients",BatchNum);
  Kokkos::DualView<Scalar*[6][6],ExeSpace> stiffness("stiffness",BatchNum);
  Kokkos::DualView<Scalar*[6], ExeSpace> stress("stress",BatchNum);

  auto deformation_gradient_h = deformation_gradient.h_view;
  for(int i=0; i<BatchNum; ++i)
    for(int ei=0; ei<3; ++ei)
      for(int ej=0; ej<3; ++ej)
        deformation_gradient_h(i,ei,ej) = cases[0].pd.deformationGradient[ei*3+ej];
  deformation_gradient.modify<Kokkos::HostSpace>();

  Kokkos::Timer timer;
  double time = timer.seconds();
  timer.reset();
  bool result = batched_analysis.run(deformation_gradient, stress);
  stress.sync<Kokkos::HostSpace>();
  auto stress_h = stress.h_view;

  //for(int i=0; i<6; ++i)
  //{
  //  std::cout<<stress_h(0, i)<<" ";
  //}
  std::cout<<std::endl;
  apf::Matrix3x3 strss;
  stressToMat(0, stress.h_view, strss);
  std::cout<<strss<<std::endl;
  double time2 = timer.seconds();
  std::cout<<"Took: "<<time2 << " seconds."<<std::endl;
  amsi::freeAnalysis();
  return 0;
}
