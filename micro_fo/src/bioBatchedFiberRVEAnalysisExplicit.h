#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H__
#include <memory>
#include <vector>
#include "bioBatchedRVEAnalysis.h"
#include "bioFiberNetwork.h"
#include "bioMicroFOParams.h"
#include "bioMicroTypeDefinitions.h"
//#include "bioViewOfView.h"
#include "bioPackedData.h"
#include "bioBatchedFiberRVEAnalysisExplicit_impl.h"
#include <apf.h>
#include <apfMesh.h> // for count owned
#include <apfMeshIterator.h> // iterator for RVE Functions
namespace bio
{
  template <typename Scalar = bio::Scalar,
            typename LocalOrdinal = bio::LocalOrdinal,
            typename ExeSpace = Kokkos::DefaultExecutionSpace>
  class BatchedFiberRVEAnalysisExplicit : public BatchedRVEAnalysis<Scalar>
  {
    private:
    using LO = LocalOrdinal;
    // the packed types correspond to data which is packed to include multiple
    // RVEs. Each in the packed data corresponds to the data for that RVE
    using PackedScalarType = PackedData<Scalar, LocalOrdinal, ExeSpace>;
    using PackedOrdinalType = PackedData<LocalOrdinal, LocalOrdinal, ExeSpace>;
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

    int countBoundaryDof(apf::Mesh* mesh)
    {
      auto bgn = amsi::apfMeshIterator(network, 0);
      decltype(bgn) end = amsi::apfEndIterator(mesh);
        getBoundaryVerts(this->rve.get(),
                         this->mFiberNetwork->getNetworkMesh(),
                         bgn,
                         end,
                         RVE::side::all,
                         std::back_inserter(this->bnd_nds[sd]));
    }



    public:
    //class TagBatchLoopSingle{};
    //class TagBatchLoopTeamOuterWhile {};
    // class TagBatchLoopTeamInnerWhile {};
    BatchedFiberRVEAnalysisExplicit(
        std::vector<std::unique_ptr<FiberNetwork>> fiber_networks,
        std::vector<std::unique_ptr<MicroSolutionStrategy>>
            solution_strategies)
    {
      // Get the size of each analysis' componenents to setup packed data
      std::vector<LO> dof_counts;
      std::vector<LO> fixed_dof_counts;
      std::vector<LO> element_counts;
      // we know that each packed data will have the same number of rows
      // as the number of fiber networks so reserve space for each row
      dof_counts.reserve(fiber_networks.size());
      element_counts.reserve(fiber_networks.size());
      fixed_dof_counts.reserve(fiber_networks.size());
      for(auto & fiber_network : fiber_networks)
      {
        dof_counts.push_back(fiber_network->getDofCount());
        // since we are dealing with trusses all elements are just lines
        element_counts.push_back(apf::countOwned(fiber_network->getNetworkMesh(),1));
        // 3 DOF per fixed node (ux,uy,uz). Note that the number of fixed nodes in each
        // analysis should not change every time it is run (only the values change), so
        // we do not need to update the sizes each time we run the analsis
        fixed_dof_counts.push_back(3*bnd_nds[RVE::all].size());

      }
      // initialize analysis PackedData

    };

    virtual bool run(const std::vector<DeformationGradient> & dfmGrds,
                     std::vector<Scalar[6]> sigma,
                     bool update_coords = true) final
    {
      return false;
    };
    virtual void computeMaterialStiffness(std::vector<Scalar[36]> C) final {};
  };
  // the actual explicit instantionations can be found in the associated cc file
  extern template class BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal>;
}  // namespace bio
#endif
