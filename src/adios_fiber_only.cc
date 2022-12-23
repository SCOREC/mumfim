#include <apf.h>
#include <mpi.h>
#include <iostream>
#include "mumfim/microscale/MicroFOParams.h"
#include "mumfim/microscale/FiberNetworkLibrary.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <mumfim/microscale/BatchedFiberRVEAnalysisExplicit.h>
#include <memory>
#include <PCU.h>
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
  amsi::MPI mpi{argc, argv, MPI_COMM_WORLD};
#ifdef MICRO_USING_PETSC
  las::initPETScLAS(&argc, &argv, MPI_COMM_WORLD);
#endif
  Kokkos::ScopeGuard kokkos(argc, argv);
  if(argc != 3)
  {
    std::cerr<<"Usage: "<<argv[0]<<" job.yaml num_rves"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  PCU_Switch_Comm(MPI_COMM_SELF);
  auto BatchNum = std::atoi(argv[2]);
  std::vector<mumfim::MicroCase> cases;
  mumfim::loadMicroFOFromYamlFile(argv[1], cases);
  mumfim::printMicroFOCase(cases[0]);
  std::string file_name = cases[0].pd.meshFile;
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  mumfim::FiberNetworkLibrary network_library;
  auto fiber_network =
      network_library.load(file_name, file_name + ".params", 0, 0);

  std::vector<std::shared_ptr<const mumfim::FiberNetwork>> fiber_networks;
  std::vector<std::shared_ptr<const mumfim::MicroSolutionStrategy>>
      solution_strategies;
  std::shared_ptr<mumfim::MicroSolutionStrategy> shared_case{std::move(cases[0].ss)};
  fiber_networks.reserve(BatchNum);
  for(int i=0; i<BatchNum; ++i)
  {
    fiber_networks.push_back(fiber_network);
    solution_strategies.push_back(shared_case);
  }
  
  using ExeSpace = Kokkos::DefaultExecutionSpace;
  // using ExeSpace = Kokkos::Serial;
  // using Scalar = float;//mumfim::Scalar;
  using Scalar = mumfim::Scalar;
  using Ordinal = mumfim::LocalOrdinal;
  mumfim::BatchedFiberRVEAnalysisExplicit<Scalar, Ordinal, ExeSpace>
      batched_analysis(std::move(fiber_networks),
                       std::move(solution_strategies));
  Kokkos::DualView<Scalar*[3][3],ExeSpace> deformation_gradient("deformation gradients",BatchNum);
  Kokkos::DualView<Scalar*[6][6],ExeSpace> stiffness("stiffness",BatchNum);
  Kokkos::DualView<Scalar*[6], ExeSpace> stress("stress",BatchNum);
  Kokkos::DualView<Scalar * [3][3], ExeSpace> orientation_tensor(
      "orientation tensor", BatchNum);
  Kokkos::DualView<Scalar * [3], ExeSpace> normal("normal", BatchNum);
  auto deformation_gradient_h = deformation_gradient.h_view;
  for(int i=0; i<BatchNum; ++i)
  {
    normal.h_view(i, 0) = 1;
    normal.h_view(i, 1) = 0;
    normal.h_view(i, 2) = 0;
    for(int ei=0; ei<3; ++ei)
    {
      for(int ej=0; ej<3; ++ej)
      {
        deformation_gradient_h(i,ei,ej) = cases[0].pd.deformationGradient[ei*3+ej];
      }
    }
  }
  normal.modify<Kokkos::HostSpace>();
  deformation_gradient.modify<Kokkos::HostSpace>();
  Kokkos::Timer timer;
  double ornt_time1 = timer.seconds();
  batched_analysis.compute3DOrientationTensor(orientation_tensor);
  orientation_tensor.sync<Kokkos::HostSpace>();
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      printf("%f ", orientation_tensor.h_view(0, i, j));
    }
    printf("\n");
  }
  printf("\n");
  batched_analysis.compute2DOrientationTensor(normal, orientation_tensor);
  orientation_tensor.sync<Kokkos::HostSpace>();
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      printf("%f ", orientation_tensor.h_view(0, i, j));
    }
    printf("\n");
  }
  printf("\n");
  double ornt_time2 = timer.seconds();
  double time1 = timer.seconds();
  auto success = !batched_analysis.run(deformation_gradient, stress);
  stress.sync<Kokkos::HostSpace>();
  auto stress_h = stress.h_view;
  std::cout<<std::endl;
  apf::Matrix3x3 strss;
  stressToMat(0, stress_h, strss);
  std::cout<<strss<<std::endl;
  batched_analysis.computeMaterialStiffness(stiffness);
  stiffness.sync<Kokkos::HostSpace>();
  auto stiffness_h = stiffness.h_view;
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      std::cout << stiffness_h(0, i, j) << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
  double time2 = timer.seconds();
  std::cout << "Took: " << time2 - time1 << " seconds." << std::endl;
  std::cout << "Orientation Computation Took: " << ornt_time2 - ornt_time1
            << " seconds." << std::endl;
  return success;
}
