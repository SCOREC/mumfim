#ifndef MUMFIM_BATCHED_RVE_ANALYSIS_H
#define MUMFIM_BATCHED_RVE_ANALYSIS_H
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <array>
#include <vector>
#include "MicroFOParams.h"
#include "mumfim/microscale/MicroTypeDefinitions.h"
#include "mumfim/microscale/ContinuumMehcanics.h"
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
    virtual ~BatchedRVEAnalysis()= default;;
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
    // this should be dSigma/dE, where Sigma is the cauchy stress, and E is the
    // PK2 stress. This should have 36 components due to the symmetry in Sigma
    // and E
    // default implementation uses finite difference method
    virtual void computeMaterialStiffness(
        Kokkos::DualView<Scalar * [6][6], ExeSpace> C) {
      using HostMemorySpace = typename decltype(C)::host_mirror_space;
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
