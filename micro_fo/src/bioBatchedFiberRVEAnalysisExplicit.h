#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H__
#include <memory>
#include <vector>
#include "bioBatchedRVEAnalysis.h"
#include "bioFiberNetwork.h"
#include "bioMicroFOParams.h"
#include "bioMicroTypeDefinitions.h"
#include "bioViewOfView.h"
namespace bio
{
  template <typename Scalar = bio::Scalar,
            typename LocalOrdinal = bio::LocalOrdinal,
            typename ExeSpace = Kokkos::DefaultExecutionSpace>
  class BatchedFiberRVEAnalysisExplicit : public BatchedRVEAnalysis<Scalar>
  {
    private:
    using LO = LocalOrdinal;
    // Simulation state vectors
    ViewOfView<LO, ExeSpace> connectivity;
    ViewOfView<Scalar, ExeSpace> original_coordinates;
    ViewOfView<Scalar, ExeSpace> current_coordinates;
    ViewOfView<Scalar, ExeSpace> displacement;
    ViewOfView<Scalar, ExeSpace> velocity;
    ViewOfView<Scalar, ExeSpace> acceleration;
    ViewOfView<Scalar, ExeSpace> force_total;
    ViewOfView<Scalar, ExeSpace> force_internal;
    ViewOfView<Scalar, ExeSpace> force_external;
    ViewOfView<Scalar, ExeSpace> force_damping;
    ViewOfView<Scalar, ExeSpace> nodal_mass;
    ViewOfView<Scalar, ExeSpace> original_length;
    ViewOfView<Scalar, ExeSpace> current_length;
    // Only deal with fiber networks with uniform material properties
    ViewOfView<Scalar, ExeSpace> fiber_elastic_modulus;
    ViewOfView<Scalar, ExeSpace> fiber_area;
    ViewOfView<Scalar, ExeSpace> fiber_density;
    ViewOfView<Scalar, ExeSpace> viscous_damping_coefficient;

    // Solver parameters
    ViewOfView<Scalar, ExeSpace> critical_time_scale_factor;
    ViewOfView<Scalar, ExeSpace> energy_check_epsilon;

    // Boundary Condition data
    ViewOfView<LO, ExeSpace> displacement_boundary_dof;
    ViewOfView<Scalar, ExeSpace> displacement_boundary_values;



    public:
    virtual bool run(const std::vector<DeformationGradient> & dfmGrds,
                     std::vector<Scalar[6]> sigma,
                     bool update_coords = true) final;
    virtual void computeMaterialStiffness(std::vector<Scalar[36]> C) final;
    BatchedFiberRVEAnalysisExplicit(
        std::vector<std::unique_ptr<FiberNetwork>> fiber_networks,
        std::vector<std::unique_ptr<MicroSolutionStrategy>>
            solution_strategies);
  };
  // the actual explicit instantionations can be found in the associated cc file
  extern template class BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal>;
}  // namespace bio
// implementation file for templated free functions
#include "bioBatchedFiberRVEAnalysisExplicit_impl.h"
namespace bio
{
  // Constructor needs to initialize all of the Views
  template <typename Scalar, typename LocalOrdinal, typename ExeSpace>
  BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal, ExeSpace>::
      BatchedFiberRVEAnalysisExplicit(
          std::vector<std::unique_ptr<FiberNetwork>> fiber_networks,
          std::vector<std::unique_ptr<MicroSolutionStrategy>>
              solution_strategies)
  {
    
  }
  // templated function definitions
  template <typename Scalar, typename LocalOrdinal, typename ExeSpace>
  bool BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal, ExeSpace>::run(
      const std::vector<DeformationGradient> & dfmGrds,
      std::vector<Scalar[6]> sigma,
      bool update_coords)
  {
    auto result = impl::BatchedExplicit<impl::BatchLoopType::SINGLE, LO, Scalar,
                          ExeSpace>::run(connectivity, original_coordinates,
                                         current_coordinates, displacement,
                                         velocity, acceleration, force_total,
                                         force_internal, force_external,
                                         force_damping, nodal_mass,
                                         original_length, current_length,
                                         fiber_elastic_modulus, fiber_area,
                                         fiber_density, viscous_damping_coefficient,
                                         critical_time_scale_factor, energy_check_epsilon,
                                         displacement_boundary_dof, displacement_boundary_values);
    return result;
  };
  template <typename Scalar, typename LocalOrdinal, typename ExeSpace>
  void BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal, ExeSpace>::
      computeMaterialStiffness(std::vector<Scalar[36]> C){};
}  // namespace bio
#endif
