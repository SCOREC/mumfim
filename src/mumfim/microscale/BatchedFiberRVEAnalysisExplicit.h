#ifndef MUMFIM_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H
#define MUMFIM_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H
#include <memory>
#include <vector>
#include "BatchedRVEAnalysis.h"
#include "FiberNetwork.h"
#include "MicroFOParams.h"
#include "MicroTypeDefinitions.h"
#include "PackedData.h"
//#include "bioBatchedFiberRVEAnalysisExplicitSerialOuterLoop.h"
#include <apf.h>
#include <apfMesh.h>          // for count owned
#include <apfMeshIterator.h>  // iterator for RVE Functions
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include "mumfim/microscale/batched_analysis/BatchedComputeOrientationTensor.h"
#include "mumfim/microscale/batched_analysis/BatchedFiberRVEAnalysisExplicitTeamOuterLoop.h"
#include "mumfim/microscale/batched_analysis/BatchedMesh.h"
#include "mumfim/microscale/batched_analysis/BatchedReorderMesh.h"
#include "RVE.h"
namespace mumfim
{
  // void updateRVECoords(RVE &rve, const DeformationGradient &
  // incremental_deformation_gradient);
  template <typename Scalar = mumfim::Scalar,
            typename LocalOrdinal = mumfim::LocalOrdinal,
            typename ExeSpace = Kokkos::DefaultExecutionSpace>
  class BatchedFiberRVEAnalysisExplicit
      : public BatchedRVEAnalysis<Scalar, LocalOrdinal, ExeSpace>
  {
    private:
    using LO = LocalOrdinal;
    // the packed types correspond to data which is packed to include multiple
    // RVEs. Each in the packed data corresponds to the data for that RVE
    using ConnectivityType = PackedData<LocalOrdinal*[2], ExeSpace>;
    using DofDataType = PackedData<Scalar*[3], ExeSpace>;
    using VertDataType = PackedData<Scalar*, ExeSpace>;
    using OrdinalVertDataType = PackedData<LocalOrdinal*, ExeSpace>;
    using ElementDataType = PackedData<Scalar*, ExeSpace>;

    using PackedScalarType = PackedData<Scalar*, ExeSpace>;
    using HostMemorySpace = typename PackedScalarType::host_mirror_space;
    using DeviceMemorySpace = typename PackedScalarType::memory_space;
    using DofCopyType = typename DofDataType::DeviceViewType;
    // helper class typedefs
    using MeshFunctionType =
        BatchedApfMeshFunctions<Scalar, LocalOrdinal, HostMemorySpace>;
    // the layout of all of the arrays should be the same
    using ViewLayout = typename ConnectivityType::traits::array_layout;
    using OrientationTensorType =
        OrientationTensor<Scalar, LocalOrdinal, DofDataType, ConnectivityType, ExeSpace>;
    // Simulation state vectors
    ConnectivityType connectivity;
    DofDataType original_coordinates;
    DofDataType current_coordinates;
    DofDataType displacement;
    DofDataType velocity;
    DofDataType force_total;
    DofDataType force_internal;
    VertDataType nodal_mass;
    ElementDataType original_length;
    ElementDataType current_length;
    // Only deal with fiber networks with uniform material properties
    PackedScalarType fiber_elastic_modulus;
    PackedScalarType fiber_area;
    PackedScalarType fiber_density;
    PackedScalarType viscous_damping_coefficient;
    // Solver parameters
    PackedScalarType critical_time_scale_factor;
    // PackedScalarType energy_check_epsilon;
    // Boundary Condition data
    OrdinalVertDataType displacement_boundary_vert;
    DofCopyType displacement_copy;
    Kokkos::DualView<Scalar *, ExeSpace> current_volume;
    Kokkos::DualView<Scalar *, ExeSpace> scale_factor;
    std::vector<RVE> rves;
    OrientationTensorType orientation_tensor;

    public:
    BatchedFiberRVEAnalysisExplicit(
        std::vector<std::shared_ptr<const FiberNetwork>> fiber_networks,
        std::vector<std::shared_ptr<const MicroSolutionStrategy>>
            solution_strategies) : BatchedRVEAnalysis<Scalar,LocalOrdinal,ExeSpace>(fiber_networks.size())
    {
      if (fiber_networks.size() != solution_strategies.size())
      {
        std::cerr << "There must be a solution strategy for each fiber network!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      using RowCountType = typename PackedScalarType::IndexViewType;
      // Get the size of each analysis' componenents to setup packed data
      RowCountType dof_counts("dof counts", fiber_networks.size());
      RowCountType fixed_dof_counts("fixed dof counts", fiber_networks.size());
      RowCountType element_counts("element counts", fiber_networks.size());
      RowCountType material_counts("material counts", fiber_networks.size());
      RowCountType connectivity_counts("material counts",
                                       fiber_networks.size());
      // get the host data views
      auto dof_counts_h = dof_counts.h_view;
      auto fixed_dof_counts_h = fixed_dof_counts.h_view;
      auto element_counts_h = element_counts.h_view;
      auto material_counts_h = material_counts.h_view;
      auto connectivity_counts_h = connectivity_counts.h_view;
      // mark the host side as modified
      dof_counts.template modify<HostMemorySpace>();
      fixed_dof_counts.template modify<HostMemorySpace>();
      element_counts.template modify<HostMemorySpace>();
      material_counts.template modify<HostMemorySpace>();
      connectivity_counts.template modify<HostMemorySpace>();
      // fill the row counts data
      std::vector<std::vector<apf::MeshEntity *>> boundary_verts;
      rves.reserve(fiber_networks.size());
      boundary_verts.reserve(boundary_verts.size());
      current_volume = Kokkos::DualView<Scalar *, ExeSpace>(
          "current_volume", fiber_networks.size());
      deep_copy(current_volume.template view<ExeSpace>(), 1.0);
      current_volume.template modify<ExeSpace>();
      scale_factor = Kokkos::DualView<Scalar *, ExeSpace>(
          "scale_factor", fiber_networks.size());
      scale_factor.template modify<HostMemorySpace>();
      for (size_t i = 0; i < fiber_networks.size(); ++i)
      {
        scale_factor.h_view(i) = fiber_networks[i]->getScaleConversion();
        dof_counts_h(i) = fiber_networks[i]->getDofCount()/3;
        // since we are dealing with trusses all elements are just lines
        element_counts_h(i) =
            apf::countOwned(fiber_networks[i]->getNetworkMesh(), 1);
        // we assume that we are using trusses, so there are 2 vertices for each
        // element
        connectivity_counts_h(i) = element_counts_h(i);
        // the standard RVE we use is 0.5x0.5x0.5
        rves.push_back(RVE(0.5));
        boundary_verts.push_back({});
        // get a list of all of the boundary vertices
        auto bgn =
            amsi::apfMeshIterator(fiber_networks[i]->getNetworkMesh(), 0);
        decltype(bgn) end =
            amsi::apfEndIterator(fiber_networks[i]->getNetworkMesh());
        getBoundaryVerts(&(rves.back()), fiber_networks[i]->getNetworkMesh(),
                         bgn, end, RVE::side::all,
                         std::back_inserter(boundary_verts.back()));
        // 3 DOF per fixed node (ux,uy,uz). Note that the number of fixed nodes
        // in each analysis should not change every time it is run (only the
        // values change), so we do not need to update the sizes each time we
        // run the analsis
        fixed_dof_counts_h(i) = boundary_verts.back().size();
        // fixed_dof_counts_h(i) = boundary_verts.back().size();
        // currently we only deal with one material property for all fibers
        // eventually this will be changed
        material_counts_h(i) = 1;
      }
      // we just reuse the material counts since both the residual
      // and the material properties only have one value per network for now
      // eventually, the residual will need its own count because there
      // will be one material per element
      auto residual_counts = material_counts;
      // initialize analysis PackedData
      // here we are currently assuming that we are using truss elements,
      // so the number of vertices on each element is 2
      connectivity = ConnectivityType(connectivity_counts);
      original_coordinates = DofDataType(dof_counts);
      current_coordinates = DofDataType(dof_counts);
      displacement = DofDataType(dof_counts);
      displacement_copy = DofCopyType(displacement.template getAllRows<ExeSpace>());
      displacement_copy = Kokkos::create_mirror(
          ExeSpace(), displacement.template getAllRows<ExeSpace>());
      velocity = DofDataType(dof_counts);
      force_total = DofDataType(dof_counts);
      force_internal = DofDataType(dof_counts);
      nodal_mass = VertDataType(dof_counts);
      original_length = ElementDataType(element_counts);
      current_length = ElementDataType(element_counts);
      // currently there is only one material property,
      fiber_elastic_modulus = PackedScalarType(material_counts);
      fiber_area = PackedScalarType(material_counts);
      fiber_density = PackedScalarType(material_counts);
      viscous_damping_coefficient = PackedScalarType(material_counts);
      critical_time_scale_factor = PackedScalarType(residual_counts);
      // Boundary Condition data
      displacement_boundary_vert = OrdinalVertDataType(fixed_dof_counts);
      // fill arrays that will be constant
      // we have to do this loop in the host space because most of the
      // data, and the mesh is not designed to go on the GPU.
      for (size_t i = 0; i < fiber_networks.size(); ++i)
      {
        auto connectivity_row =
            connectivity.template getRow<HostMemorySpace>(i);
        MeshFunctionType::getConnectivity(fiber_networks[i]->getNetworkMesh(),
                                          connectivity_row, 1);
        auto original_coordinates_row =
            original_coordinates.template getRow<HostMemorySpace>(i);
        MeshFunctionType::getFieldValues(
            fiber_networks[i]->getNetworkMesh(),
            fiber_networks[i]->getNetworkMesh()->getCoordinateField(),
            original_coordinates_row);
        auto displacement_boundary_vert_row =
            displacement_boundary_vert.template getRow<HostMemorySpace>(i);
        apf::NaiveOrder(fiber_networks[i]->getUNumbering());
        MeshFunctionType::getFixedVert(fiber_networks[i]->getUNumbering(),
                                       displacement_boundary_vert_row, boundary_verts[i]);
        auto fiber_elastic_modulus_row =
            fiber_elastic_modulus.template getRow<HostMemorySpace>(i);
        fiber_elastic_modulus_row(0) =
            fiber_networks[i]->getFiberReaction(0).getYoungModulus();
        auto fiber_area_row = fiber_area.template getRow<HostMemorySpace>(i);
        fiber_area_row(0) =
            fiber_networks[i]->getFiberReaction(0).getFiberArea();
        auto fiber_density_row =
            fiber_density.template getRow<HostMemorySpace>(i);
        fiber_density_row(0) =
            fiber_networks[i]->getFiberReaction(0).getFiberDensity();
        auto viscous_damping_coefficient_row =
            viscous_damping_coefficient.template getRow<HostMemorySpace>(i);
        viscous_damping_coefficient_row(0) =
            static_cast<const MicroSolutionStrategyExplicit *>(
                solution_strategies[i].get())
                ->visc_damp_coeff;
        auto critical_time_scale_factor_row =
            critical_time_scale_factor.template getRow<HostMemorySpace>(i);
        critical_time_scale_factor_row(0) =
            static_cast<const MicroSolutionStrategyExplicit *>(
                solution_strategies[i].get())
                ->crit_time_scale_factor;
        auto nodal_mass_row = nodal_mass.template getRow<HostMemorySpace>(i);
        MeshFunctionType::getNodalMass(fiber_networks[i]->getNetworkMesh(),
                                       fiber_density_row(0), fiber_area_row(0),
                                       nodal_mass_row);
        // reorder the mesh
        ReorderMesh<Scalar, LocalOrdinal, ViewLayout, Kokkos::Serial> reorder(
                 nodal_mass_row.extent(0), displacement_boundary_vert_row);
        reorder.createPermutationArray();
        reorder.applyPermutationToConnectivity(connectivity_row);
        reorder.applyPermutationToCoordinates(original_coordinates_row);
        reorder.applyPermutationToNodalValue(nodal_mass_row);
        reorder.applyPermutationToFixedVert();
      }
      // mark the packed data as modified
      connectivity.template modify<HostMemorySpace>();
      original_coordinates.template modify<HostMemorySpace>();
      nodal_mass.template modify<HostMemorySpace>();
      fiber_elastic_modulus.template modify<HostMemorySpace>();
      fiber_area.template modify<HostMemorySpace>();
      fiber_density.template modify<HostMemorySpace>();
      viscous_damping_coefficient.template modify<HostMemorySpace>();
      critical_time_scale_factor.template modify<HostMemorySpace>();
      displacement_boundary_vert.template modify<HostMemorySpace>();
      // we need to initialize the current coordinates for the orientaiton
      // tensor since we want to compute this before we do any analysis runs
      Kokkos::deep_copy(
          current_coordinates.template getAllRows<HostMemorySpace>(),
          original_coordinates.template getAllRows<HostMemorySpace>());
      current_coordinates.template modify<HostMemorySpace>();
      orientation_tensor =
          OrientationTensorType(connectivity, current_coordinates, TEAM_SIZE);
    }
    virtual bool run(
        Kokkos::DualView<Scalar * [3][3], ExeSpace> deformation_gradients,
        Kokkos::DualView<Scalar * [6], ExeSpace> sigma,
        bool update_coords = true) final
    {
      displacement.template sync<ExeSpace>();
      displacement.template modify<ExeSpace>();
      if (!update_coords)
      {
        Kokkos::deep_copy(displacement_copy,
                          displacement.template getAllRows<ExeSpace>());
      }
      else
      {
        updateVolume<ExeSpace>(deformation_gradients, current_volume);
      }
      // apply the incremental deformation gradient to the boundaries and the
      // boundary_vert_values
      applyIncrementalDeformationToDisplacement<ExeSpace>(
          rves.size(), deformation_gradients, current_coordinates,
          displacement);
      // sync all data to the execution space
      connectivity.template sync<ExeSpace>();
      original_coordinates.template sync<ExeSpace>();
      current_coordinates.template sync<ExeSpace>();
      displacement.template sync<ExeSpace>();
      velocity.template sync<ExeSpace>();
      force_total.template sync<ExeSpace>();
      force_internal.template sync<ExeSpace>();
      nodal_mass.template sync<ExeSpace>();
      original_length.template sync<ExeSpace>();
      current_length.template sync<ExeSpace>();
      fiber_elastic_modulus.template sync<ExeSpace>();
      fiber_area.template sync<ExeSpace>();
      fiber_density.template sync<ExeSpace>();
      viscous_damping_coefficient.template sync<ExeSpace>();
      critical_time_scale_factor.template sync<ExeSpace>();
      displacement_boundary_vert.template sync<ExeSpace>();

      Kokkos::deep_copy(velocity.template getAllRows<ExeSpace>(),0);
      // Run the specific implementation we are interested in
      auto result = TeamOuterLoop<Scalar, LO, ExeSpace>::run(rves.size(),
          connectivity, original_coordinates, current_coordinates,
          displacement, velocity, force_total, force_internal,
          nodal_mass, original_length,
          current_length, fiber_elastic_modulus, fiber_area,
          fiber_density, viscous_damping_coefficient,
          critical_time_scale_factor, displacement_boundary_vert, TEAM_SIZE);
      // compute the stress from the force and displacement vectors
      computeCauchyStress<ExeSpace>(displacement_boundary_vert,
                                    current_coordinates, force_internal, sigma,
                                    current_volume, scale_factor);
      if (!update_coords)
      {
        Kokkos::deep_copy(displacement.template getAllRows<ExeSpace>(),
                          displacement_copy);
      }
      else
      {
        Kokkos::deep_copy(this->current_stress_.template view<ExeSpace>(), sigma.template view<ExeSpace>());
        this->current_stress_.template modify<ExeSpace>();
      }
      return result;
    }
    virtual void computeMaterialStiffness(
        Kokkos::DualView<Scalar * [6][6], ExeSpace> C) final
    {
      if (this->current_stress_.extent(0) != C.extent(0))
      {
        std::cerr << "stiffness matrix must have the number of rves as the "
                     "first dimension"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      constexpr Scalar h = 1E-5;
      this->current_stress_.template sync<ExeSpace>();
      auto current_stress_d = this->current_stress_.template view<ExeSpace>();
      Kokkos::DualView<Scalar * [6], ExeSpace> sigma1("C", C.extent(0));
      auto sigma1_d = sigma1.template view<ExeSpace>();
      Kokkos::DualView<Scalar * [3][3], ExeSpace> Fappd("Fappd", 1);
      auto Fappd_h = Fappd.template view<HostMemorySpace>();
      // auto Fappd_d = Fappd.template view<ExeSpace>();
      auto C_d = C.template view<ExeSpace>();
      Scalar D1 = sqrt(1.0 / (1 - 2 * h));
      constexpr Scalar l1 = 1;
      constexpr Scalar l2 = (1 + h) / (1 - h * h);
      constexpr Scalar l3 = (1 - h) / (1 - h * h);
      Scalar l2pl3 = 0.5 * (sqrt(l2) + sqrt(l3));
      Scalar l2ml3 = 0.5 * (sqrt(l2) - sqrt(l3));
      // compute V from V = sqrt((I-2e)^-1)
      // V has been computed by hand
      for (int i = 0; i < 6; ++i)
      {
        Kokkos::deep_copy(Fappd_h, 0);
        switch (i)
        {
          case 0:
            Fappd_h(0, 0, 0) = D1;
            Fappd_h(0, 1, 1) = 1;
            Fappd_h(0, 2, 2) = 1;
            break;
          case 1:
            Fappd_h(0, 0, 0) = 1;
            Fappd_h(0, 1, 1) = D1;
            Fappd_h(0, 2, 2) = 1;
            break;
          case 2:
            Fappd_h(0, 0, 0) = 1;
            Fappd_h(0, 1, 1) = 1;
            Fappd_h(0, 2, 2) = D1;
            break;
          case 3:
            Fappd_h(0, 0, 0) = l1;
            Fappd_h(0, 1, 1) = l2pl3;
            Fappd_h(0, 2, 2) = l2pl3;
            Fappd_h(0, 1, 2) = l2ml3;
            Fappd_h(0, 2, 1) = l2ml3;
            break;
          case 4:
            Fappd_h(0, 0, 0) = l2pl3;
            Fappd_h(0, 1, 1) = l1;
            Fappd_h(0, 2, 2) = l2pl3;
            Fappd_h(0, 0, 2) = l2ml3;
            Fappd_h(0, 2, 0) = l2ml3;
            break;
          case 5:
            Fappd_h(0, 0, 0) = l2pl3;
            Fappd_h(0, 1, 1) = l2pl3;
            Fappd_h(0, 2, 2) = l1;
            Fappd_h(0, 0, 1) = l2ml3;
            Fappd_h(0, 1, 0) = l2ml3;
            break;
        }
        Fappd.template modify<HostMemorySpace>();
        run(Fappd, sigma1, false);
        sigma1.template sync<ExeSpace>();
        // j is row, i is column
        Kokkos::parallel_for(
            Kokkos::RangePolicy<ExeSpace>(0, current_stress_d.extent(0)),
            KOKKOS_LAMBDA(const int rve) {
              for (int j = 0; j < 6; ++j)
              {
                C_d(rve, j, i) =
                    (sigma1_d(rve, j) - current_stress_d(rve, j)) / h;
              }
            });
      }
      C.template modify<ExeSpace>();
    }
    void compute3DOrientationTensor(
        Kokkos::DualView<Scalar * [3][3], ExeSpace> omega) final
    {
      orientation_tensor.compute3D(omega);
    }
    void compute2DOrientationTensor(
        Kokkos::DualView<Scalar * [3], ExeSpace> normal,
        Kokkos::DualView<Scalar * [3][3], ExeSpace> omega) final
    {
      orientation_tensor.compute2D(normal, omega);
    }
  };
}  // namespace mumfim
#endif
