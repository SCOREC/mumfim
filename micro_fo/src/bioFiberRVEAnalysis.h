#ifndef BIO_FIBER_RVE_ANALYSIS_H_
#define BIO_FIBER_RVE_ANALYSIS_H_
#include "bioFiberNetwork.h"
#include "bioRVE.h"
#include "bioTrussIntegrator.h"
#include "bioMultiscaleRVE.h"
#include "lasSparskit.h"
#include <amsiAnalysis.h>
#include <string>
namespace bio
{
  struct FiberRVEAnalysis
  {
    FiberNetwork * fn;
    // move this out of here! requires coupling data structures as currently written hmmm...
    MultiscaleRVE * multi;
    RVE * rve;
    std::vector<apf::MeshEntity*> bnd_nds[RVE::side::all+1];
    apf::Integrator * es;
    las::Mat * k;
    las::Vec * u;
    las::Vec * f;
    las::Vec * f0;
    las::LasOps * ops;
    las::LasSolve * slv;
    apf::DynamicMatrix dx_fn_dx_rve;
    bool dx_fn_dx_rve_set;
  };
  FiberRVEAnalysis * makeFiberRVEAnalysis(FiberNetwork *, las::CSR *, las::SparskitBuffers * b = NULL);
  void destroyAnalysis(FiberRVEAnalysis *);
  class FiberRVEIteration : public amsi::Iteration
  {
  protected:
    FiberRVEAnalysis * an;
  public:
    FiberRVEIteration(FiberRVEAnalysis * a);
    void iterate();
  };
   /**
    * Fix the boundary dofs and set the
    *  rows in the mat/vec corresponing
    *  to the dof numbering to identity rows
    */
  template <typename I>
  void applyRVEBC(I bnd_bgn,
                  I bnd_end,
                  apf::Numbering * nm,
                  las::LasOps * ops,
                  las::Mat * k,
                  las::Vec * f);
  /**
   * Free all boundary dofs so they can
   *   be numbered prior to being fixed again.
   */
  template <typename I>
  void freeeRVEBC(I bnd_bgn,
                  I bnd_end,
                  apf::Numbering * num);
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma);
}
#include "bioFiberRVEAnalysis_impl.h"
#endif
