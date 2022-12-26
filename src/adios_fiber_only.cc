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
#include <adios2.h>
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

auto ReadAdiosDeformations(const char* inputfile, adios2::ADIOS& ad, MPI_Comm comm=MPI_COMM_WORLD) {
  int rank, size;
  MPI_Comm_rank(comm, &rank); 
  MPI_Comm_size(comm, &size);

  adios2::IO bpReaderIO = ad.DeclareIO("input");
  if (!bpReaderIO.InConfigFile())
 	{
 	//		// if not defined by user, we can change the default settings
 	//		// BPFile is the default engine
 			bpReaderIO.SetEngine("BP");
 			  //bpReaderIO.SetEngine("BP5");
 	//		//bpReaderIO.SetParameters({{"num_threads", "2"}});
 	//		bpReaderIO.SetParameter("OpenAsFile", "true");
 
 	//		// ISO-POSIX file is the default transport
 	//		// Passing parameters to the transport
 	//		bpReaderIO.AddTransport("File", {{"verbose", "4"}});
 	}
	adios2::Engine bpReader = bpReaderIO.Open(inputfile, adios2::Mode::Read, comm);
	
	
	auto v = bpReaderIO.InquireVariable<double>("F");
	auto shape = v.Shape();
	
	size_t N = shape[0]/size;
	adios2::Dims start{rank*N, 0};
	//fixme count is 
	auto cnt = (rank==(size-1))?shape[0]-rank*N:N;
	adios2::Dims count{cnt, shape[1]};
  Kokkos::View<double*[3][3], Kokkos::LayoutRight, Kokkos::HostSpace> data("F", cnt);
  v.SetSelection({start,count});
  bpReader.Get(v, data.data(), adios2::Mode::Sync);
	bpReader.Close();
  return std::pair {shape[0], data};
}

int main(int argc, char * argv[])
{

  amsi::MPI mpi{argc, argv, MPI_COMM_WORLD};
#ifdef MICRO_USING_PETSC
  las::initPETScLAS(&argc, &argv, MPI_COMM_WORLD);
#endif
  Kokkos::ScopeGuard kokkos(argc, argv);
  adios2::ADIOS ad(MPI_COMM_WORLD);
  if(argc != 4)
  {
    std::cerr<<"Usage: "<<argv[0]<<" deformation.bp output.bp network_name.txt"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  const char *inputfile = argv[1];
  const char *outputfile = argv[2];
  const char *network_name = argv[3];
  //const char *num_realizations = std::atoi(argv[4]);

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  auto [num_total, deformation_gradients_in] = ReadAdiosDeformations(inputfile, ad, MPI_COMM_WORLD);
  size_t N = num_total/size;
	auto BatchNum = (rank==(size-1))?num_total-rank*N:N;


  auto solution_strategy = std::make_shared<mumfim::MicroSolutionStrategyExplicit>();
  solution_strategy->cnvgTolerance=1E-6;
  solution_strategy->slvrTolerance=1E-6;
  solution_strategy->slvrType=mumfim::SolverType::Explicit;
  solution_strategy->oscPrms={};
  solution_strategy->print_history_frequency=0;
  solution_strategy->print_field_frequency=0;
  solution_strategy->print_field_by_num_frames=false;
  solution_strategy->visc_damp_coeff = 0.7;
  solution_strategy->load_time=5;
  solution_strategy->total_time=100;
  solution_strategy->serial_gpu_cutoff=0;
  solution_strategy->crit_time_scale_factor=0.9;
  solution_strategy->energy_check_eps=0.01;
  solution_strategy->ampType = mumfim::AmplitudeType::SmoothStep;
  
  std::vector<std::shared_ptr<const mumfim::MicroSolutionStrategy>>
      solution_strategies(BatchNum, solution_strategy);

  // this is needed for loading the libraries...due to pumi
  PCU_Switch_Comm(MPI_COMM_SELF);
  mumfim::FiberNetworkLibrary network_library;
  
  auto fiber_network =
      network_library.load(network_name, std::string(network_name) + ".params", 0, 0);

  std::vector<std::shared_ptr<const mumfim::FiberNetwork>> fiber_networks(BatchNum, fiber_network);
  
  using ExeSpace = Kokkos::DefaultExecutionSpace;
  using Scalar = mumfim::Scalar;
  using Ordinal = mumfim::LocalOrdinal;
  mumfim::BatchedFiberRVEAnalysisExplicit<Scalar, Ordinal, ExeSpace>
      batched_analysis(std::move(fiber_networks),
                       std::move(solution_strategies));
  Kokkos::DualView<Scalar*[3][3],ExeSpace> deformation_gradient("deformation gradients",BatchNum);
  Kokkos::DualView<Scalar*[6][6],ExeSpace> stiffness("stiffness",BatchNum);
  Kokkos::DualView<Scalar*[6], ExeSpace> stress("stress",BatchNum);
  //Kokkos::DualView<Scalar * [3][3], ExeSpace> orientation_tensor("orientation tensor", BatchNum);


  auto deformation_gradient_h = deformation_gradient.h_view;
  for(int i=0; i<BatchNum; ++i)
  {
    for(int ei=0; ei<3; ++ei)
    {
      for(int ej=0; ej<3; ++ej)
      {
        deformation_gradient_h(i,ei,ej) = deformation_gradients_in(i, ei, ej);
      }
    }
  }
  deformation_gradient.modify<Kokkos::HostSpace>();

  //deformation_gradient.modify<Kokkos::HostSpace>();
  //batched_analysis.compute3DOrientationTensor(orientation_tensor);
  //orientation_tensor.sync<Kokkos::HostSpace>();
  auto success = !batched_analysis.run(deformation_gradient, stress);
  stress.sync<Kokkos::HostSpace>();
  auto stress_h = stress.h_view;
  batched_analysis.computeMaterialStiffness(stiffness);
  stiffness.sync<Kokkos::HostSpace>();
  auto stiffness_h = stiffness.h_view;

  PCU_Switch_Comm(MPI_COMM_WORLD);

  adios2::IO bpIO = ad.DeclareIO("output");
  if (!bpIO.InConfigFile())
 	{
 			bpIO.SetEngine("BP");
 	}
  //auto orientation_tensor_var = bpIO.DefineVariable<Scalar>("orientation_tensor", {num_total,3,3}, {rank*N,0,0}, {BatchNum, 3,3}, true);
  auto stress_var = bpIO.DefineVariable<Scalar>("stresss", {num_total,6}, {rank*N,0}, {BatchNum, 6}, true);
  auto stiffness_var = bpIO.DefineVariable<Scalar>("stiffness", {num_total,6,6}, {rank*N,0,0}, {BatchNum, 6,6}, true);
	
	adios2::Engine bpWriter = bpIO.Open(outputfile, adios2::Mode::Write, MPI_COMM_WORLD);
  //bpWriter.Put(orientation_tensor_var, orientation_tensor.h_view.data());
  bpWriter.Put(stress_var, stress_h.data());
  bpWriter.Put(stiffness_var, stiffness_h.data());

	bpWriter.Close();
  return success;
}
