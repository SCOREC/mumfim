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
#include <numeric>
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
 			bpReaderIO.SetEngine("BP");
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
  if(argc != 5)
  {
    std::cerr<<"Usage: "<<argv[0]<<" deformation.bp output.bp network_name.txt write_field_data"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  const char *inputfile = argv[1];
  const char *outputfile = argv[2];
  const char *network_name = argv[3];
  bool write_field_data = std::atoi(argv[4]);
  std::cout<<"Writing field data: "<<std::boolalpha<<write_field_data<<"\n";
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
  
  std::cout<<"Initializing analysis\n";
  std::vector<std::shared_ptr<const mumfim::MicroSolutionStrategy>>
      solution_strategies(BatchNum, solution_strategy);

  // this is needed for loading the libraries...due to pumi
  PCU_Switch_Comm(MPI_COMM_SELF);
  mumfim::FiberNetworkLibrary network_library;
  
  auto fiber_network =
      network_library.load(network_name, std::string(network_name) + ".params", 0, 0);

  std::vector<std::shared_ptr<const mumfim::FiberNetwork>> fiber_networks(BatchNum, fiber_network);
  std::cout <<"Networks loaded\n";
  
  using ExeSpace = Kokkos::DefaultExecutionSpace;
  using Scalar = mumfim::Scalar;
  using Ordinal = mumfim::LocalOrdinal;
  std::cout <<"Creating adios2 IO\n";
  adios2::IO bpIO = ad.DeclareIO("output");
  std::cout <<"Creating adios2 Engine\n";
  if (!bpIO.InConfigFile())
 	{
 			bpIO.SetEngine("BP");
 	}
  std::cout <<"Open adios2 Engine\n";
	adios2::Engine bpWriter = bpIO.Open(outputfile, adios2::Mode::Write, MPI_COMM_WORLD);
  std::cout<<"Create Batched analysis\n";
  bpWriter.BeginStep();


  mumfim::BatchedFiberRVEAnalysisExplicit<Scalar, Ordinal, ExeSpace>
      batched_analysis(std::move(fiber_networks),
                       std::move(solution_strategies));
  Kokkos::DualView<Scalar*[3][3],ExeSpace> deformation_gradient("deformation gradients",BatchNum);
  Kokkos::DualView<Scalar*[6][6],ExeSpace> stiffness("stiffness",BatchNum);
  Kokkos::DualView<Scalar*[6], ExeSpace> stress("stress",BatchNum);
  Kokkos::DualView<Scalar * [3][3], ExeSpace> orientation_tensor("orientation tensor", BatchNum);

  std::cout<<"Load deformatin gradient\n";
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

  std::cout<<"Running\n";
  auto success = !batched_analysis.run(deformation_gradient, stress);
  stress.sync<Kokkos::HostSpace>();
  std::cout<<"Done running\n";
  if(write_field_data) {
    auto connectivity_packed = batched_analysis.GetConnectivityArray();
    connectivity_packed.sync<Kokkos::HostSpace>();
    auto coordinates_packed = batched_analysis.GetCoordinates();
    coordinates_packed.sync<Kokkos::HostSpace>();
    auto displacement_packed = batched_analysis.GetDisplacement();
    displacement_packed.sync<Kokkos::HostSpace>();
    auto velocity_packed = batched_analysis.GetVelocity();
    velocity_packed.sync<Kokkos::HostSpace>();
    auto force_packed = batched_analysis.GetForce();
    force_packed.sync<Kokkos::HostSpace>();
    std::cout <<"Writing field data\n";
    uint64_t num_elem = batched_analysis.GetConnectivityArray().getAllRows<ExeSpace>().extent(0);
    std::cout<<"Num elements "<<num_elem<<"\n";
    uint64_t num_verts = batched_analysis.GetCoordinates().getAllRows<ExeSpace>().extent(0);
    std::cout<<"Num verts "<<num_verts<<"\n";
    // copy data from device to host no because stiffness computation will change data
    Kokkos::View<uint64_t*[2], Kokkos::LayoutRight, Kokkos::HostSpace> connectivity_h("connectivity", num_elem);
    Kokkos::deep_copy(connectivity_h,connectivity_packed.getAllRows<Kokkos::HostSpace>());
    Kokkos::View<double*[3], Kokkos::LayoutRight, Kokkos::HostSpace> coordinates_h("coordinates", num_verts);
    Kokkos::deep_copy(coordinates_h, coordinates_packed.getAllRows<Kokkos::HostSpace>());
    Kokkos::View<double*[3], Kokkos::LayoutRight, Kokkos::HostSpace> displacement_h("displacement", num_verts);
    Kokkos::deep_copy(displacement_h, displacement_packed.getAllRows<Kokkos::HostSpace>());
    Kokkos::View<double*[3], Kokkos::LayoutRight, Kokkos::HostSpace> velocity_h("velocity", num_verts);
    Kokkos::deep_copy(velocity_h, velocity_packed.getAllRows<Kokkos::HostSpace>());
    Kokkos::View<double*[3], Kokkos::LayoutRight, Kokkos::HostSpace> force_h("force", num_verts);
    Kokkos::deep_copy(force_h, force_packed.getAllRows<Kokkos::HostSpace>());
    std::cout <<"Copy done field data\n";

    auto num_rows = connectivity_packed.getNumRows();
    assert(num_rows == BatchNum);

    std::cout <<"Getting Offsets\n";
    std::vector<uint64_t> vert_offsets(BatchNum+1);
    std::vector<uint64_t> elem_offsets(BatchNum+1);
    // need to get connectivity vert_offsets
    // One offset per network corresponds to each row
    vert_offsets[0] = 0;
    elem_offsets[0] = 0;
    for(int i=1; i<BatchNum+1; ++i) {
      vert_offsets[i] = coordinates_packed.getRow<Kokkos::HostSpace>(i-1).extent(0);
      elem_offsets[i] = connectivity_packed.getRow<Kokkos::HostSpace>(i-1).extent(0);
    }
    std::cout <<"local scan\n";
    std::inclusive_scan(std::next(vert_offsets.begin()), vert_offsets.end(), std::next(vert_offsets.begin()));
    std::inclusive_scan(std::next(elem_offsets.begin()), elem_offsets.end(), std::next(elem_offsets.begin()));
    assert(vert_offsets.back() == num_verts);
    assert(elem_offsets.back() == num_elem);
    std::cout <<"MPI scan\n";
    uint64_t global_vert_offset =0, global_elem_offset=0;
    MPI_Exscan(&num_verts, &global_vert_offset, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Exscan(&num_elem, &global_elem_offset, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    printf("Rank: %d, Globl Vert: %ld, Global Elem: %ld\n", rank, global_vert_offset, global_elem_offset);
    // iterate over the connectivity data and load it into the adios connectivity
    Kokkos::View<uint64_t*[3], Kokkos::LayoutRight, Kokkos::HostSpace> adios_connectivity("adios connectivity", num_elem);
    assert(vert_offsets.size() == elem_offsets.size());
    for(int i=0; i<elem_offsets.size()-1; ++i) {
      auto start = elem_offsets[i];
      auto end = elem_offsets[i+1];
      std::cout<<start<<" "<<end<<" "<<num_elem<<"\n";
      assert(end <= num_elem);
      auto offset = vert_offsets[i];
      // need to make the offsets unique for each entry in the connectivity array
      for(int j=start; j<end; ++j) {
        adios_connectivity(j, 0) = 2;
        // don't need global vert offset because adios2 reader requires use of local arrays
        adios_connectivity(j,1) = connectivity_h(j,0)+offset;
        adios_connectivity(j,2) = connectivity_h(j,1)+offset;
      }
    }
    std::cout << "Connectivity filled\n";
    // Specify which network each vertex is associated with. This will be used in paraview to filter data
    // corresponding to the specific network
    Kokkos::View<uint64_t*, Kokkos::LayoutRight, Kokkos::HostSpace> network_id_h("network id", num_verts);
    for(int i=0; i<vert_offsets.size()-1; ++i) {
      auto start = vert_offsets[i];
      auto end = vert_offsets[i+1];
      for(int j=start; j<end; ++j) {
        network_id_h(j) = rank*N+i;
      }
    }
    std::cout<<"MPI Allreduce\n";
    uint64_t num_total_verts,num_total_elem;
    MPI_Allreduce(&num_verts, &num_total_verts,1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&num_elem, &num_total_elem,1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    std::cout<<"define adios variables\n";
    //auto connectivity_var = bpIO.DefineVariable<uint64_t>("connectivity", {num_total_elem, 3},{global_elem_offset, 0},{num_elem, 3},true);
    //auto vertices_var = bpIO.DefineVariable<Scalar>("vertices", {num_total_verts,3}, {global_vert_offset,0}, {num_verts, 3}, true);
    //auto displacement_var = bpIO.DefineVariable<Scalar>("displacement", {num_total_verts,3}, {global_vert_offset,0}, {num_verts, 3}, true);
    //auto velocity_var = bpIO.DefineVariable<Scalar>("velocity", {num_total_verts,3}, {global_vert_offset,0}, {num_verts, 3}, true);
    //auto force_var = bpIO.DefineVariable<Scalar>("force", {num_total_verts,3}, {global_vert_offset,0}, {num_verts, 3}, true);
    //auto network_id_var = bpIO.DefineVariable<uint64_t>("network", {num_total_verts}, {global_vert_offset}, {num_verts}, true);
    
    // need to use local arrays
    auto connectivity_var = bpIO.DefineVariable<uint64_t>("connectivity", {},{},{num_elem, 3},true);
    auto vertices_var = bpIO.DefineVariable<Scalar>("vertices", {}, {}, {num_verts, 3}, true);
    auto displacement_var = bpIO.DefineVariable<Scalar>("displacement", {}, {}, {num_verts, 3}, true);
    auto velocity_var = bpIO.DefineVariable<Scalar>("velocity", {}, {}, {num_verts, 3}, true);
    auto force_var = bpIO.DefineVariable<Scalar>("force", {}, {}, {num_verts, 3}, true);
    auto network_id_var = bpIO.DefineVariable<uint64_t>("network", {}, {}, {num_verts}, true);
    std::cout << "Putting data\n"; 

    bpWriter.Put(connectivity_var, adios_connectivity.data());
    bpWriter.Put(vertices_var, coordinates_h.data());
    bpWriter.Put(displacement_var, displacement_h.data());
    bpWriter.Put(velocity_var, velocity_h.data());
    bpWriter.Put(force_var, force_h.data());
    bpWriter.Put(network_id_var, network_id_h.data());
    std::cout<<"Put Type\n";
    // vtk line is type 3
    constexpr int32_t  line_type = 3;

    auto types_var = bpIO.DefineVariable<int32_t>("types");
    bpWriter.Put(types_var, line_type);
    std::cout<<"Putting attibutes\n";
    const std::string unstructureGridSchema = R"(                                                                                                                                                          
          <VTKFile type="UnstructuredGrid">
            <UnstructuredGrid>
              <Piece>
                <Points>
                  <DataArray Name="vertices" />
                </Points>
                <Cells>
                  <DataArray Name="connectivity" />
                  <DataArray Name="types" />
                </Cells>
                <PointData>
                  <DataArray Name="displacement" />
                  <DataArray Name="velocity" />
                  <DataArray Name="force" />
                  <DataArray Name="network" />
                </PointData>
              </Piece>
            </UnstructuredGrid>
          </VTKFile>)";

    bpIO.DefineAttribute("vtk.xml", unstructureGridSchema);
    
    std::cout<<"Done doing field data\n";
    bpWriter.PerformPuts();
  }


  batched_analysis.compute3DOrientationTensor(orientation_tensor);
  orientation_tensor.sync<Kokkos::HostSpace>();
  

  auto stress_h = stress.h_view;
  batched_analysis.computeMaterialStiffness(stiffness);
  stiffness.sync<Kokkos::HostSpace>();
  auto stiffness_h = stiffness.h_view;

  PCU_Switch_Comm(MPI_COMM_WORLD);

  bpIO.DefineAttribute("Network", std::string(network_name));
  auto deformation_var = bpIO.DefineVariable<Scalar>("F", {num_total,3,3}, {rank*N,0,0}, {BatchNum, 3,3}, true);
  auto orientation_var = bpIO.DefineVariable<Scalar>("orientation", {num_total,3,3}, {rank*N,0,0}, {BatchNum, 3,3}, true);
  auto stress_var = bpIO.DefineVariable<Scalar>("stress", {num_total,6}, {rank*N,0}, {BatchNum, 6}, true);
  auto stiffness_var = bpIO.DefineVariable<Scalar>("stiffness", {num_total,6,6},{rank*N,0,0}, {BatchNum, 6,6}, true);
  auto volume_var = bpIO.DefineVariable<Scalar>("volume", {num_total},{rank*N}, {BatchNum}, true);
  auto strain_energy_var = bpIO.DefineVariable<Scalar>("strain_energy", {num_total},{rank*N}, {BatchNum}, true);
  // Note that the data is stored in layout_left for GPU, so we must copy it into layout right to write it out properly for adios2
  Kokkos::View<double*[3][3], Kokkos::LayoutRight, Kokkos::HostSpace> right_deformation_gradient_h("Fh", BatchNum);
  Kokkos::deep_copy(right_deformation_gradient_h, deformation_gradient_h);

  Kokkos::View<double*[3][3], Kokkos::LayoutRight, Kokkos::HostSpace> right_orientation_h("orientation", BatchNum);
  Kokkos::deep_copy(right_orientation_h, orientation_tensor.h_view);

  Kokkos::View<double*[6], Kokkos::LayoutRight, Kokkos::HostSpace> right_stress_h("stress", BatchNum);
  Kokkos::deep_copy(right_stress_h, stress_h);

  Kokkos::View<double*[6][6], Kokkos::LayoutRight, Kokkos::HostSpace> right_stiffness_h("stiffness", BatchNum);
  Kokkos::deep_copy(right_stiffness_h, stiffness_h);

  Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> volume_h("volume", BatchNum);
  Kokkos::deep_copy(volume_h, batched_analysis.GetVolume().d_view);
  Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> strain_energy_h("strain energy", BatchNum);
  Kokkos::deep_copy(strain_energy_h, batched_analysis.GetStrainEnergy());
  std::cerr<<"Putting small vars\n";
  bpWriter.Put(deformation_var, right_deformation_gradient_h.data());
  bpWriter.Put(stress_var, right_stress_h.data());
  bpWriter.Put(orientation_var, right_orientation_h.data());
  bpWriter.Put(stiffness_var, right_stiffness_h.data());
  bpWriter.Put(volume_var, volume_h.data());
  bpWriter.Put(strain_energy_var, strain_energy_h.data());
  
  std::cerr<<"Closing\n";
  bpWriter.EndStep();
	bpWriter.Close();
  return success;
}
