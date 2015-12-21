#ifndef MICRO_FO_ANALYSIS_H_
#define MICRO_FO_ANALYSIS_H_
#include "numNewtonRaphson.h"
#include "FiberNetwork.h"
#include "RVE.h"
#include "TrussIntegrator.h"
#include <lasSparskit.h>
#include <string>
namespace bio
{
  struct FiberRVEAnalysis
  {
    FiberNetwork * fn;
    RVE * rve;
    apf::Integrator * itgr;
    apf::FieldOp * writesol;
    apf::FieldOp * accumsol;
    las::skMat * k;
    las::skVec u;
    las::skVec f;
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
}
#endif
