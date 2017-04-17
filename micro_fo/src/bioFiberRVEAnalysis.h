#ifndef BIO_FIBER_RVE_ANALYSIS_H_
#define BIO_FIBER_RVE_ANALYSIS_H_
#include "bioFiberNetwork.h"
#include "bioRVE2.h"
#include "bioTrussIntegrator.h"
#include "lasSparskit.h"
#include <amsiAnalysis.h>
#include <string>
namespace bio
{
  struct FiberRVEAnalysis
  {
    FiberNetwork * fn;
    RVE * rve;
    amsi::ElementalSystem * es
    las::Mat * k;
    las::Vec u;
    las::Vec f;
  };
  FiberRVEAnalysis * makeAnalysis(const std::string & fnm);
  void destroyAnalysis(FiberRVEAnalysis *);
  class FiberRVEIteration : public num::Iteration
  {
  protected:
    FiberRVEAnalysis * an;
  public:
    FiberRVEIteration(FiberRVEAnalysis * a);
    void iterate();
  };
  class FiberRVEConvergence : public num::Convergence
  {
  protected:
    FiberRVEAnalysis * an;
    double eps;
    double resid_im;
  public:
    FiberRVEConvergence(FiberRVEAnalysis * a, double e = 1e-8);
    bool converged();
    double & epsilon() {return eps;}
  };
  double calcStiffness(FiberRVEAnalysis * fra);
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3 & sigma);
  void tdYdXr(FiberRVEAnalysis * fra);
  void calcFEMJacob(FiberRVEAnalysis * fra);
}
#endif
