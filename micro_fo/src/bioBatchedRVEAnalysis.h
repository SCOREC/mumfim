#ifndef BIO_BATCHED_RVE_ANALYSIS_H__
#define BIO_BATCHED_RVE_ANALYSIS_H__
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <array>
#include <vector>
#include "bioMicroFOParams.h"
namespace bio
{
  template <typename Scalar, typename Ordinal, typename ExeSpace>
  class BatchedRVEAnalysis
  {
    protected:
    // TODO convert the data representation to a kokkos view
    Kokkos::DualView<Scalar * [6], ExeSpace> current_stress_;
    Ordinal num_rves_;

    public:
    BatchedRVEAnalysis(Ordinal num_rves) : num_rves_(num_rves) {
      current_stress_ = Kokkos::DualView<Scalar * [6], ExeSpace>("current_stress", num_rves);
    };
    virtual ~BatchedRVEAnalysis(){};
    // for now we only allow the case where every analysis needs to update
    // coords, or every analysis doesn't it is possible this isn't the most
    // efficient choice, and we can re-evaluate this at a later time
    // virtual bool run(const std::vector<DeformationGradient> & dfmGrds,
    //                 std::vector<Scalar[6]> sigma, bool update_coords=true) =
    //                 0;
    virtual bool run(Kokkos::DualView<Scalar * [3][3], ExeSpace> dfmGrds,
                     Kokkos::DualView<Scalar * [6], ExeSpace> sigma,
                     bool update_coords = true) = 0;
    // computes the material stiffness tensor at the current deformation state
    // this should be dSigma/dE, where Sigma is the cauchy stress, and E is the
    // PK2 stress. This should have 36 components due to the symmetry in Sigma
    // and E
    virtual void computeMaterialStiffness(
        Kokkos::DualView<Scalar * [6][6], ExeSpace> C) = 0;
    virtual void compute3DOrientationTensor(
        Kokkos::DualView<Scalar * [3][3], ExeSpace> omega)
    {
      omega.template modify<Kokkos::HostSpace>();
      Kokkos::deep_copy(omega.h_view, 0);
    }
    virtual void compute2DOrientationTensor(
        Kokkos::DualView<Scalar * [3], ExeSpace> normal,
        Kokkos::DualView<Scalar * [3][3], ExeSpace> omega)
    {
      omega.template modify<Kokkos::HostSpace>();
      Kokkos::deep_copy(omega.h_view, 0);
    }
    // RVEAnalysis(const RVEAnalysis & an);
    // RVEAnalysis();
  };
}  // namespace bio
#endif
