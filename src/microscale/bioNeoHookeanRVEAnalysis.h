#ifndef BIO_NEOHOOKEAN_RVE_ANALYSIS_H_
#define BIO_NEOHOOKEAN_RVE_ANALYSIS_H_
#include "bioMicroFOParams.h"
#include "bioMultiscaleMicroFOParams.h"
#include "bioRVEAnalysis.h"
#include <apfMatrix.h>
#include <apfDynamicMatrix.h>
namespace bio {
class NeoHookeanRVEAnalysis : public RVEAnalysis
{
  protected:
    apf::Matrix3x3 F;
    apf::Matrix3x3 F_old;
    double detF;
    apf::DynamicMatrix leftCauchyGreen;
    double ShearModulus;
    double lambda;
    virtual void computeCauchyStress(double sigma[6]);
  public:
    NeoHookeanRVEAnalysis(double youngs_modulus, double poisson_ratio)
      : ShearModulus(youngs_modulus / (2.0 * (1.0 + poisson_ratio)))
      , lambda((2.0 * ShearModulus * poisson_ratio) / (1.0 - 2.0 * poisson_ratio))
      {
        F = apf::Matrix3x3(1,0,0,
                           0,1,0,
                           0,0,1);
      }
    virtual bool run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords=true) final;
    virtual void computeMaterialStiffness(double C[36]) final;
};
  std::unique_ptr<NeoHookeanRVEAnalysis> initNeoHookeanRVEAnalysisFromMultiscale(micro_fo_params & prm);
}

#endif
