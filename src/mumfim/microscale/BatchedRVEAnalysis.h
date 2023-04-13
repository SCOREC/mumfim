#ifndef MUMFIM_BATCHED_RVE_ANALYSIS_H
#define MUMFIM_BATCHED_RVE_ANALYSIS_H
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <array>
#include <vector>
#include "MicroFOParams.h"
#include "mumfim/microscale/MicroTypeDefinitions.h"
#include "mumfim/microscale/ContinuumMechanics.h"
namespace mumfim
{
  template <typename Scalar, typename Ordinal, typename ExeSpace>
  class BatchedRVEAnalysis
  {
    protected:
    // TODO convert the data representation to a kokkos view
    Kokkos::DualView<Scalar * [6], ExeSpace> current_stress_;
    Ordinal num_rves_;

    public:
    using memory_space = typename ExeSpace::memory_space;
    using exe_space = ExeSpace;

    // TODO get rid of DualView in the interface
    explicit BatchedRVEAnalysis(Ordinal num_rves) : num_rves_(num_rves) {
      current_stress_ = Kokkos::DualView<Scalar * [6], ExeSpace>("current_stress", num_rves);
    };
    virtual ~BatchedRVEAnalysis()= default;
    // for now we only allow the case where every analysis needs to update
    // coords, or every analysis doesn't it is possible this isn't the most
    // efficient choice, and we can re-evaluate this at a later time
    // virtual bool run(const std::vector<DeformationGradient> & dfmGrds,
    //                 std::vector<Scalar[6]> sigma, bool update_coords=true) =
    //                 0;
    virtual bool run(Kokkos::DualView<Scalar * [3][3], ExeSpace> dfmGrds,
                     Kokkos::DualView<Scalar * [6], ExeSpace> sigma,
                     bool update_coords = true) = 0;
    /* Accept the state from the last run*/
    virtual void accept() {}
    // computes the material stiffness tensor at the current deformation state
    // this can be computed by F@F@D@F@F where F is deformation tensor and 
    // D= dPK2/dE
    virtual void computeMaterialStiffness(
        Kokkos::DualView<Scalar * [6][6], ExeSpace> C) = 0;
    virtual void compute3DOrientationTensor(
        Kokkos::DualView<Scalar * [3][3], ExeSpace> omega)
    {
      omega.template modify<Kokkos::HostSpace>();
      Kokkos::deep_copy(omega.h_view, 0);
    }
    virtual void compute2DOrientationTensor(
        Kokkos::DualView<Scalar * [3], ExeSpace> /*normal*/,
        Kokkos::DualView<Scalar * [3][3], ExeSpace> omega)
    {
      omega.template modify<Kokkos::HostSpace>();
      Kokkos::deep_copy(omega.h_view, 0);
    }
    Ordinal GetNumRVEs() const { return this->num_rves_; }
    // RVEAnalysis(const RVEAnalysis & an);
    // RVEAnalysis();
  };

  template <typename Analysis>
  struct BatchedAnalysisGetPK2StressFunc
  {
    using memory_space = typename Analysis::memory_space;
    using exe_space = typename Analysis::exe_space;

    explicit BatchedAnalysisGetPK2StressFunc(Analysis & analysis)
        : analysis_(analysis)
        , F_increment_("deformation grad increment", analysis_.GetNumRVEs())
        , stress_("stress", analysis_.GetNumRVEs())
    {
    }

    // That way we can work with codes that use either incrementalor total
    // approach
    auto operator()(Kokkos::View<Scalar * [3][3], memory_space> F,
                    Kokkos::View<Scalar * [3][3], memory_space> F_increment,
                    bool update_coords = false) noexcept
    -> Kokkos::View<Scalar * [6], memory_space>
    {
      assert(F_increment.extent(0) == analysis_.GetNumRVEs());
      assert(F.extent(0) == analysis_.GetNumRVEs());
      Kokkos::deep_copy(F_increment_.d_view, F_increment);
      F_increment_.template modify<exe_space>();
      analysis_.run(F_increment_, stress_, update_coords);
      stress_.template sync<exe_space>();
      // Batched analysis computes the cauchy stress, we need to convert this to
      // the PK2 stress
      ConvertCauchyToPK2<exe_space>(F, stress_.d_view);
      return stress_.d_view;
    }

    Analysis & analysis_;
    Kokkos::DualView<Scalar * [3][3], memory_space> F_increment_;
    Kokkos::DualView<Scalar * [6], memory_space> stress_;
  };

}  // namespace mumfim
#endif
