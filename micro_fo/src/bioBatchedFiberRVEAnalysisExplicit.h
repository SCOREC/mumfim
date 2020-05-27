#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H__
#include <memory>
#include <vector>
#include "bioBatchedRVEAnalysis.h"
#include "bioFiberNetwork.h"
#include "bioMicroFOParams.h"
#include "bioMicroTypeDefinitions.h"
#include "bioPackedData.h"
#include "bioBatchedFiberRVEAnalysisExplicit_impl.h"
#include <apf.h>
#include <apfMesh.h> // for count owned
#include <apfMeshIterator.h> // iterator for RVE Functions
#include "bioRVE.h"
namespace bio
{
  void updateRVECoords(RVE &rve, const DeformationGradient & incremental_deformation_gradient);

  template <typename Scalar = bio::Scalar,
            typename LocalOrdinal = bio::LocalOrdinal,
            typename ExeSpace = Kokkos::DefaultExecutionSpace>
  class BatchedFiberRVEAnalysisExplicit : public BatchedRVEAnalysis<Scalar,ExeSpace>
  {
    private:
    using LO = LocalOrdinal;
    // the packed types correspond to data which is packed to include multiple
    // RVEs. Each in the packed data corresponds to the data for that RVE
    using PackedScalarType = PackedData<Scalar, LocalOrdinal, ExeSpace>;
    using PackedOrdinalType = PackedData<LocalOrdinal, LocalOrdinal, ExeSpace>;
    using HostMemorySpace = typename PackedScalarType::host_mirror_space;
    using DeviceMemorySpace = typename PackedScalarType::memory_space;
    // Simulation state vectors
    PackedOrdinalType connectivity;
    PackedScalarType original_coordinates;
    PackedScalarType current_coordinates;
    PackedScalarType displacement;
    PackedScalarType velocity;
    PackedScalarType acceleration;
    PackedScalarType force_total;
    PackedScalarType force_internal;
    PackedScalarType force_external;
    PackedScalarType force_damping;
    PackedScalarType nodal_mass;
    PackedScalarType original_length;
    PackedScalarType current_length;
    PackedScalarType residual; // solutoin residual
    // Only deal with fiber networks with uniform material properties
    PackedScalarType fiber_elastic_modulus;
    PackedScalarType fiber_area;
    PackedScalarType fiber_density;
    PackedScalarType viscous_damping_coefficient;

    // Solver parameters
    PackedScalarType critical_time_scale_factor;
    // PackedScalarType energy_check_epsilon;

    // Boundary Condition data
    PackedOrdinalType displacement_boundary_dof;
    PackedScalarType displacement_boundary_values;

    std::vector<RVE> rves;


    public:
    //class TagBatchLoopSingle{};
    //class TagBatchLoopTeamOuterWhile {};
    // class TagBatchLoopTeamInnerWhile {};
    BatchedFiberRVEAnalysisExplicit(
        std::vector<std::unique_ptr<FiberNetwork>> fiber_networks,
        std::vector<std::unique_ptr<MicroSolutionStrategy>>
            solution_strategies)
    {
      if(fiber_networks.size() != solution_strategies.size())
      {
        std::cerr<<"There must be a solution strategy for each fiber network!"<<std::endl;
        std::exit(EXIT_FAILURE);
      }
      using RowCountType = typename PackedScalarType::IndexViewType;
      // Get the size of each analysis' componenents to setup packed data
      RowCountType dof_counts("dof counts", fiber_networks.size());
      RowCountType fixed_dof_counts("fixed dof counts", fiber_networks.size());
      RowCountType element_counts("element counts", fiber_networks.size());
      RowCountType material_counts("material counts", fiber_networks.size());
      RowCountType connectivity_counts("material counts", fiber_networks.size());
      RowCountType mass_counts("mass counts", fiber_networks.size());
      // get the host data views
      auto dof_counts_h = dof_counts.h_view;
      auto fixed_dof_counts_h = fixed_dof_counts.h_view;
      auto element_counts_h = element_counts.h_view;
      auto material_counts_h = material_counts.h_view;
      auto connectivity_counts_h = connectivity_counts.h_view;
      auto mass_counts_h = mass_counts.h_view;
      // mark the host side as modified
      dof_counts.template modify<HostMemorySpace>();
      fixed_dof_counts.template modify<HostMemorySpace>();
      element_counts.template modify<HostMemorySpace>();
      material_counts.template modify<HostMemorySpace>();
      connectivity_counts.template modify<HostMemorySpace>();
      mass_counts.template modify<HostMemorySpace>();
      // fill the row counts data
      std::vector<std::vector<apf::MeshEntity*>> boundary_verts;
      rves.reserve(fiber_networks.size());
      boundary_verts.reserve(boundary_verts.size());
      for(size_t i=0; i<fiber_networks.size(); ++i)
      {
        dof_counts_h(i)=fiber_networks[i]->getDofCount();
        mass_counts_h(i)=dof_counts_h(i)/3;
        // since we are dealing with trusses all elements are just lines
        element_counts_h(i) = apf::countOwned(fiber_networks[i]->getNetworkMesh(),1);
        // we assume that we are using trusses, so there are 2 vertices for each element
        connectivity_counts_h(i) = 2*element_counts_h(i);
        // the standard RVE we use is 0.5x0.5x0.5
        rves.push_back(RVE(0.5));
        boundary_verts.push_back({});
        // get a list of all of the boundary vertices
        auto bgn = amsi::apfMeshIterator(fiber_networks[i]->getNetworkMesh(), 0);
        decltype(bgn) end = amsi::apfEndIterator(fiber_networks[i]->getNetworkMesh());
        getBoundaryVerts(&(rves.back()),
                         fiber_networks[i]->getNetworkMesh(),
                         bgn,
                         end,
                         RVE::side::all,
                         std::back_inserter(boundary_verts.back()));
        // 3 DOF per fixed node (ux,uy,uz). Note that the number of fixed nodes in each
        // analysis should not change every time it is run (only the values change), so
        // we do not need to update the sizes each time we run the analsis
        fixed_dof_counts_h(i) = 3*boundary_verts.back().size();
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
      connectivity = PackedOrdinalType(connectivity_counts);
      original_coordinates = PackedScalarType(dof_counts);
      current_coordinates = PackedScalarType(dof_counts);
      displacement = PackedScalarType(dof_counts);
      velocity = PackedScalarType(dof_counts);
      acceleration = PackedScalarType(dof_counts);
      force_total = PackedScalarType(dof_counts);
      force_internal = PackedScalarType(dof_counts);
      force_external = PackedScalarType(dof_counts);
      force_damping = PackedScalarType(dof_counts);
      nodal_mass = PackedScalarType(mass_counts);
      original_length = PackedScalarType(element_counts);
      current_length = PackedScalarType(element_counts);
      // currently there is only one material property,
      residual = PackedScalarType(residual_counts);
      fiber_elastic_modulus = PackedScalarType(material_counts);
      fiber_area = PackedScalarType(material_counts);
      fiber_density = PackedScalarType(material_counts);
      viscous_damping_coefficient = PackedScalarType(material_counts);
      critical_time_scale_factor = PackedScalarType(residual_counts);
      // Boundary Condition data
      displacement_boundary_dof = PackedOrdinalType(fixed_dof_counts);
      displacement_boundary_values = PackedScalarType(fixed_dof_counts);
      // fill arrays that will be constant 
      // we have to do this loop in the host space because most of the
      // data, and the mesh is not designed to go on the GPU.
      for(size_t i=0; i<fiber_networks.size(); ++i)
      {
        auto connectivity_row = connectivity.template getRow<HostMemorySpace>(i);
        getConnectivity(fiber_networks[i]->getNetworkMesh(), connectivity_row, 1);
        auto original_coordinates_row = original_coordinates.template getRow<HostMemorySpace>(i);
        getFieldValues(fiber_networks[i]->getNetworkMesh(),
                       fiber_networks[i]->getNetworkMesh()->getCoordinateField(),
                       original_coordinates_row);
        auto displacement_boundary_dof_row = displacement_boundary_dof.template getRow<HostMemorySpace>(i);
        apf::NaiveOrder(fiber_networks[i]->getUNumbering());
        getFixedDof(fiber_networks[i]->getUNumbering(), displacement_boundary_dof_row, boundary_verts[i]);
        auto fiber_elastic_modulus_row = fiber_elastic_modulus.template getRow<HostMemorySpace>(i);
        fiber_elastic_modulus_row(0) = fiber_networks[i]->getFiberReaction(0).E;
        auto fiber_area_row = fiber_area.template getRow<HostMemorySpace>(i);
        fiber_area_row(0) = fiber_networks[i]->getFiberReaction(0).fiber_area;
        auto fiber_density_row = fiber_density.template getRow<HostMemorySpace>(i);
        fiber_density_row(0) = fiber_networks[i]->getFiberReaction(0).fiber_density;
        auto viscous_damping_coefficient_row = viscous_damping_coefficient.template getRow<HostMemorySpace>(i);
        viscous_damping_coefficient_row(0) = static_cast<MicroSolutionStrategyExplicit*>(solution_strategies[i].get())->visc_damp_coeff;
        auto critical_time_scale_factor_row = critical_time_scale_factor.template getRow<HostMemorySpace>(i);
        critical_time_scale_factor_row(0) = static_cast<MicroSolutionStrategyExplicit*>(solution_strategies[i].get())->crit_time_scale_factor;

        auto nodal_mass_row = nodal_mass.template getRow<HostMemorySpace>(i);
        getNodalMass(fiber_networks[i]->getNetworkMesh(), fiber_density_row(0),
                     fiber_area_row(0), nodal_mass_row);
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
      displacement_boundary_dof.template modify<HostMemorySpace>();
    };
    bool runSingleLoop()
    {
      // sync all data to the execution space
      connectivity.template sync<ExeSpace>();
      original_coordinates.template sync<ExeSpace>();
      current_coordinates.template sync<ExeSpace>();
      displacement.template sync<ExeSpace>();
      velocity.template sync<ExeSpace>();
      acceleration.template sync<ExeSpace>();
      force_total.template sync<ExeSpace>();
      force_internal.template sync<ExeSpace>();
      force_external.template sync<ExeSpace>();
      force_damping.template sync<ExeSpace>();
      nodal_mass.template sync<ExeSpace>();
      original_length.template sync<ExeSpace>();
      current_length.template sync<ExeSpace>();
      residual.template sync<ExeSpace>();
      fiber_elastic_modulus.template sync<ExeSpace>();
      fiber_area.template sync<ExeSpace>();
      fiber_density.template sync<ExeSpace>();
      viscous_damping_coefficient.template sync<ExeSpace>();
      critical_time_scale_factor.template sync<ExeSpace>();
      displacement_boundary_dof.template sync<ExeSpace>();
      displacement_boundary_values.template sync<ExeSpace>();
      for(size_t i=0; i<rves.size(); ++i)
      {
        //// get the device views of the data
        auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
        auto original_coordinates_row = original_coordinates.template getRow<ExeSpace>(i);
        auto current_coordinates_row = current_coordinates.template getRow<ExeSpace>(i);
        auto displacement_row = displacement.template getRow<ExeSpace>(i);
        auto velocity_row = velocity.template getRow<ExeSpace>(i);
        auto acceleration_row = acceleration.template getRow<ExeSpace>(i);
        auto force_total_row = force_total.template getRow<ExeSpace>(i);
        auto force_internal_row = force_internal.template getRow<ExeSpace>(i);
        auto force_external_row = force_external.template getRow<ExeSpace>(i);
        auto force_damping_row = force_damping.template getRow<ExeSpace>(i);
        auto nodal_mass_row = nodal_mass.template getRow<ExeSpace>(i);
        auto original_length_row = original_length.template getRow<ExeSpace>(i);
        auto current_length_row = current_length.template getRow<ExeSpace>(i);
        auto residual_row = residual.template getRow<ExeSpace>(i);

        auto displacement_boundary_dof_row = displacement_boundary_dof.template getRow<ExeSpace>(i);
        auto displacement_boundary_values_row = displacement_boundary_values.template getRow<ExeSpace>(i);

        // get the scalar values
        Scalar visc_damp_coeff = viscous_damping_coefficient.template getRow<HostMemorySpace>(i)(0);
        Scalar crit_time_scale_factor = critical_time_scale_factor.template getRow<HostMemorySpace>(i)(0);
        Scalar elastic_modulus = fiber_elastic_modulus.template getRow<HostMemorySpace>(i)(0);
        Scalar area = fiber_area.template getRow<HostMemorySpace>(i)(0);
        Scalar density = fiber_density.template getRow<HostMemorySpace>(i)(0);
        // create the loop policies based on the lengths of the various arrays
        Kokkos::RangePolicy<ExeSpace> element_policy(0, current_length_row.extent(0));
        Kokkos::RangePolicy<ExeSpace> dof_policy(0, displacement_row.extent(0));
        Kokkos::RangePolicy<ExeSpace> fixed_dof_policy(0, displacement_boundary_dof_row.extent(0));
        Scalar residual = 10.0;
        Scalar total_time = 10.0;
        Scalar current_time = 10.0;
        Scalar dt_nphalf, t_npone, t_nphalf;
        long int n_step=0;
        getElementLengths(element_policy, original_coordinates_row, connectivity_row, original_length_row);
        getCurrentCoords(dof_policy, original_coordinates_row, displacement_row, current_coordinates_row);
        getElementLengths(element_policy, current_coordinates_row, connectivity_row, current_length_row);
        auto dt_crit = getForces(element_policy, dof_policy,
                            elastic_modulus, area, density,
                            original_length_row,
                            current_length_row, 
                            velocity_row, current_coordinates_row,
                            connectivity_row, nodal_mass_row, visc_damp_coeff,
                            force_internal_row, force_external_row,
                            force_damping_row, force_total_row, residual);
        updateAcceleration(dof_policy, nodal_mass_row, force_total_row, acceleration_row);

        do
        {
          dt_nphalf = crit_time_scale_factor *dt_crit;
          t_npone = current_time + dt_nphalf;
          t_nphalf = 0.5* (current_time + t_npone);
          assert(dt_nphalf != t_npone != t_nphalf);
          updateVelocity(dof_policy, acceleration_row, (t_nphalf - current_time), velocity_row);
          updateDisplacement(dof_policy, velocity_row, dt_nphalf, displacement_row);
          applyDispBC(fixed_dof_policy, 1.0, displacement_boundary_dof_row,
                      displacement_boundary_values_row, displacement_row);
          getCurrentCoords(dof_policy, original_coordinates_row, displacement_row,
                           current_coordinates_row);
          getElementLengths(element_policy, current_coordinates_row, connectivity_row,
                            current_length_row);
          dt_crit = getForces(element_policy, dof_policy,
                              elastic_modulus, area, density,
                              original_length_row,
                              current_length_row, 
                              velocity_row, current_coordinates_row,
                              connectivity_row, nodal_mass_row, visc_damp_coeff,
                              force_internal_row, force_external_row,
                              force_damping_row, force_total_row, residual);
          updateAccelVel(dof_policy, t_npone-t_nphalf, nodal_mass_row, force_total_row,
                         acceleration_row, velocity_row);
          // the derivative and second derivative amplitudes are zero
          // since we apply the displacement boundary condition all at once
          applyAccelVelBC(fixed_dof_policy,(Scalar)0,(Scalar)0,visc_damp_coeff,displacement_boundary_dof_row,
                       displacement_boundary_values_row, nodal_mass_row, acceleration_row,
                       velocity_row, force_internal_row, force_external_row, force_damping_row,
                       force_total_row);
          current_time = t_npone;
          ++n_step;
          Kokkos::fence();
        } while ((current_time < total_time) || (residual > 1E-6));
        printf("n_step %ld, %e\n", n_step, residual);
      }
      return true;
    }

    virtual bool run(Kokkos::DualView<Scalar*[3][3], ExeSpace> deformation_gradients,
                     Kokkos::DualView<Scalar*[6], ExeSpace> sigma,
                     bool update_coords=true) final
    {
      // for coordinate update can we just set the intital coordinates
      // to the initial_coordinates+displacements if we want update?
      deformation_gradients.template sync<HostMemorySpace>();
      // apply the incremental deformation to the RVE
      for(size_t i=0; i<rves.size(); ++i)
      {
        auto deformation_gradient = Kokkos::subview(deformation_gradients, i, Kokkos::ALL, Kokkos::ALL);
        updateRVECoords(rves[i], deformation_gradient.h_view);
      }
      // apply the incremental deformation gradient to the boundaries and the boundary_dof_values
      applyIncrementalDeformationToDisplacement<ExeSpace>(deformation_gradients, original_coordinates, displacement);
      applyIncrementalDeformationToBoundary<ExeSpace>(deformation_gradients, original_coordinates, displacement_boundary_dof, displacement_boundary_values);
      // solve the problem?
      auto result = runSingleLoop();
      // compute the stress from the force and displacement vectors
      computeCauchyStress<ExeSpace>(displacement_boundary_dof, current_coordinates, force_internal, sigma);
      return result;
    };
    virtual void computeMaterialStiffness(Kokkos::DualView<Scalar*[6][6], ExeSpace> C) final {};
  };
  // the actual explicit instantionations can be found in the associated cc file
  extern template class BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal>;
}  // namespace bio
#endif
