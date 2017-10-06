#ifndef BIO_FIBER_RVE_ANALYSIS_H_
#define BIO_FIBER_RVE_ANALYSIS_H_
#include "bioFiberNetwork.h"
#include "bioRVE2.h"
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
    MultiscaleRVE * multi;
    RVE * rve;
    std::vector<apf::MeshEntity*> bnd_nds;
    apf::Integrator * es;
    las::Mat * k;
    las::Vec * u;
    las::Vec * f;
    las::LasOps * ops;
    las::LasSolve * slv;
  };
  FiberRVEAnalysis * makeFiberRVEAnalysis(FiberNetwork *,las::SparskitBuffers * b = NULL);
  void destroyAnalysis(FiberRVEAnalysis *);
  class FiberRVEIteration : public amsi::Iteration
  {
  protected:
    FiberRVEAnalysis * an;
  public:
    FiberRVEIteration(FiberRVEAnalysis * a);
    void iterate();
  };
  class FiberRVEConvergence : public amsi::Convergence
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
  /**
   * Set the force vector value associated with dofs on nodes lying on the RVE
   *  boundaries to zero.
   * @param f The force vector
   * @param rve The rve which defines which nodes are on the boundary
   * @param fn The fiber network to check for boundary nodes
   */
  template <typename I>
    void applyRVEBC(I bnd_bgn, I bnd_end, apf::Numbering * nm);
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma);
}
#include "bioFiberRVEAnalysis_impl.h"
#endif
