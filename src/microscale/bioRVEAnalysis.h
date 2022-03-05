#ifndef BIO_RVE_ANALYSIS_H_
#define BIO_RVE_ANALYSIS_H_
#include "bioMicroFOParams.h"
namespace bio
{
  class RVEAnalysis
  {
    protected:
    // this state should be updated in the run functions
    double curStress[6];

    public:
    virtual ~RVEAnalysis() {};
    // Current implementations of this function update all internal
    // variables. Be wary of this...
    virtual bool run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords=true) = 0;
    // computes the cauchy stress at the current deformation state
    //virtual void computeCauchyStress(double sigma[6]) = 0;
    // computes the material stiffness tensor at the current deformation state
    // this should be dSigma/dE, where Sigma is the cauchy stress, and E is the
    // PK2 stress. This should have 36 components due to the symmetry in Sigma and E
    virtual void computeMaterialStiffness(double C[36]);
    virtual void computeAvgVolStress(double Q[3]) { Q[0] = Q[1] = Q[2] = 0; }
    RVEAnalysis(const RVEAnalysis & an);
    RVEAnalysis();
  };
} // namespace bio
#endif
