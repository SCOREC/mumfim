#ifndef BIO_BATCHED_RVE_ANALYSIS_H__
#define BIO_BATCHED_RVE_ANALYSIS_H__
#include <vector>
#include "bioMicroFOParams.h"
#include <array>
#include <Kokkos_Core.hpp>

namespace bio
{
  template <typename Scalar>
  class BatchedRVEAnalysis
  {
    protected:
    // TODO convert the data representation to a kokkos view
    Kokkos::View<Scalar*[6]> current_stress_;

    public:
    virtual ~BatchedRVEAnalysis() {};
    // for now we only allow the case where every analysis needs to update coords,
    // or every analysis doesn't it is possible this isn't the most efficient choice,
    // and we can re-evaluate this at a later time
    virtual bool run(const std::vector<DeformationGradient> & dfmGrds,
                     std::vector<Scalar[6]> sigma, bool update_coords=true) = 0;
    // computes the material stiffness tensor at the current deformation state
    // this should be dSigma/dE, where Sigma is the cauchy stress, and E is the
    // PK2 stress. This should have 36 components due to the symmetry in Sigma and E
    virtual void computeMaterialStiffness(std::vector<Scalar[36]> C);
    virtual void computeAvgVolStress(double Q[3]) { Q[0] = Q[1] = Q[2] = 0; }
    //RVEAnalysis(const RVEAnalysis & an);
    //RVEAnalysis();
  };
} // namespace bio

#endif
