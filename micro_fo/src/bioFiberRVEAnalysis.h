#ifndef BIO_FIBER_RVE_ANALYSIS_H_
#define BIO_FIBER_RVE_ANALYSIS_H_
#include "bioFiberNetwork.h"
#include "bioRVE.h"
#include "bioTrussIntegrator.h"
#include "bioMultiscaleRVE.h"
#include <lasSparskit.h>
#include <amsiNonlinearAnalysis.h>
#include <string>
namespace bio
{
  class FiberRVEAnalysis;
  // here we have a helper class that owns the las vectors and matricies.
  // This lets us share the memory between multiple FiberRVEAnalysis instances
  // (especially copies)
  class FiberRVEAnalysisVecs {
    public:
    FiberRVEAnalysisVecs(int ndofs, las::Sparsity *csr,
                         las::SparskitBuffers *b);
    ~FiberRVEAnalysisVecs();
    // give direct access to buffers to avoid too much function call overhead
    friend class FiberRVEAnalysis;

    protected:
    las::Mat *k;
    las::Vec *u;
    las::Vec *f;
    las::Solve *slv;
  };
  class FiberRVEAnalysis
  {
    protected:
    FiberNetwork * fn;
    public:
    // move this out of here! requires coupling data structures as currently written hmmm...
    MultiscaleRVE * multi;
    RVE * rve;
    std::vector<apf::MeshEntity*> bnd_nds[RVE::side::all+1];
    apf::Integrator * es;
    apf::DynamicMatrix dx_fn_dx_rve;
    FiberRVEAnalysisVecs* vecs;
    bool dx_fn_dx_rve_set;
    double solver_eps;
    double prev_itr_factor;
    int max_cut_attempt;
    int attempt_cut_factor;
    int max_itrs;
    amsi::DetectOscillationType detect_osc_type;
    // constructors
    explicit FiberRVEAnalysis(const FiberRVEAnalysis & an);
    FiberRVEAnalysis(FiberNetwork *fn, FiberRVEAnalysisVecs *vecs,
                     micro_fo_solver &slvr, micro_fo_int_solver &slvr_int);
    ~FiberRVEAnalysis();
    las::Mat *getK() const { return vecs->k; }
    las::Vec *getU() const { return vecs->u; }
    las::Vec *getF() const { return vecs->f; }
    las::Solve *getSlv() const { return vecs->slv; }
    FiberNetwork *getFn() const { return fn; }
  };
  FiberRVEAnalysis *createFiberRVEAnalysis(FiberNetwork *fn,
                                           FiberRVEAnalysisVecs *vecs,
                                           micro_fo_solver &slvr,
                                           micro_fo_int_solver &slvr_int);
  FiberRVEAnalysis *initFromMultiscale(
      FiberNetwork *fn, FiberRVEAnalysisVecs *vecs, micro_fo_header &hdr,
      micro_fo_params &prm, micro_fo_init_data &ini, micro_fo_solver &slvr,
      micro_fo_int_solver &slvr_int);
  void destroyAnalysis(FiberRVEAnalysis *);
  FiberRVEAnalysisVecs *createFiberRVEAnalysisVecs(
      int ndofs, las::Sparsity *csr, las::SparskitBuffers *bfrs = NULL);
  void destroyFiberRVEAnalysisVecs(FiberRVEAnalysisVecs *vecs);
  /*
   * perform a deep copy of the fiber rve analysis
   */
  FiberRVEAnalysis *  copyAnalysis(FiberRVEAnalysis * an);
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
