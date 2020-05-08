#ifndef BIO_FIBER_RVE_ANALYSIS_H_
#define BIO_FIBER_RVE_ANALYSIS_H_
#include <amsiNonlinearAnalysis.h>
#include <lasSparskit.h>
#include <string>
#include "bioFiberNetwork.h"
#include "bioMassIntegrator.h"
#include "bioMicroFOParams.h"
#include "bioRVE.h"
#include "bioTrussIntegrator.h"
#include "bioMultiscaleMicroFOParams.h"
#include "bioRVEAnalysis.h"
namespace bio
{
  /* \brief create the elemental system integrator
   */
  apf::Integrator * createImplicitMicroElementalSystem(
      FiberNetwork * fn,
      las::Mat * k,
      las::Vec * f);
  class FiberRVEAnalysis;
  // here we have a helper class that owns the las vectors and matricies.
  // This lets us share the memory between multiple FiberRVEAnalysis instances
  // (especially copies)
  //
  // FIXME move this class to the static analysis class since it is not needed in the general case!
  template <typename T>
  class LinearStructs
  {
    public:
    // here we are sacrificing some type safety for a
    // consistent interface if we are using sparskit
    // the buffers should be a las::SparskitBuffers pointer
    // and they are not used for petsc, so it can be a nullptr
    LinearStructs(int ndofs,
                  double solver_tol,
                  las::Sparsity * csr,
                  void * bfrs = nullptr,
                  las::Sparsity * massCsr = nullptr);
    las::Mat * getK() const { return k; }
    las::Vec * getU() const { return u; }
    las::Vec * getF() const { return f; }
    // swaps the vector pointers if zero is true the v2 vector
    // will be zero'd
    las::Solve * getSlv() const { return slv; }
    ~LinearStructs();
    // give direct access to buffers to avoid too much function call overhead
    friend class FiberRVEAnalysis;
    friend class FiberRVEAnalysisSImplicit;
    protected:
    las::Mat * k;
    las::Vec * u;
    las::Vec * f;
    las::Solve * slv;
  };
  class FiberRVEAnalysis : public RVEAnalysis
  {
    protected:
    std::unique_ptr<FiberNetwork> mFiberNetwork;
    std::unique_ptr<MicroSolutionStrategy> mSolutionStrategy;
    RVE * rve;
    virtual void computeCauchyStress(double sigma[6]);
    // FIXME make this into a unique_ptr
    LinearStructs<las::MICRO_BACKEND> * vecs;

    public:
    std::vector<apf::MeshEntity *> bnd_nds[RVE::side::all + 1];
    apf::Integrator * es;
    double solver_eps;
    double prev_itr_factor;
    int max_cut_attempt;
    int attempt_cut_factor;
    int max_itrs;
    amsi::DetectOscillationType detect_osc_type;
    FiberRVEAnalysis(const FiberRVEAnalysis & an);
    FiberRVEAnalysis(std::unique_ptr<FiberNetwork> fn,
                     std::unique_ptr<MicroSolutionStrategy> ss);
    virtual ~FiberRVEAnalysis();
    FiberNetwork * getFn() const { return mFiberNetwork.get(); }
    RVE * getRVE() const { return rve; }
    virtual bool run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords=true) override = 0;
    virtual SolverType getAnalysisType() = 0;

    // FIXME move these functions related to linear vectors to Static class since they are not needed in
    // the explicit case
    // get the global stiffness matrix
    las::Mat * getK() const { return vecs->k; }
    // get the global displacement vector
    las::Vec * getU() const { return vecs->u; }
    // get the global force vector
    las::Vec * getF() const { return vecs->f; }
    // solve is needed for computing derivs
    las::Solve * getSlv() const { return vecs->slv; }
  };
  std::unique_ptr<FiberRVEAnalysis> createFiberRVEAnalysis(
      std::unique_ptr<FiberNetwork> fiber_network,
      std::unique_ptr<MicroSolutionStrategy> solution_strategy);
  // FIXME move to bioMultiscaleRVEAnalysis
  std::unique_ptr<FiberRVEAnalysis> initFiberRVEAnalysisFromMultiscale(
      std::unique_ptr<FiberNetwork> fiber_network,
      micro_fo_header & hdr,
      micro_fo_params & prm,
      std::unique_ptr<MicroSolutionStrategy> solution_strategy);
  LinearStructs<las::MICRO_BACKEND> * createLinearStructs(
      int ndofs,
      double solver_tol,
      las::Sparsity * csr,
      void * bfrs = nullptr,
      las::Sparsity * massSprs = nullptr);
  void destroyFiberRVEAnalysisLinearStructs(
      LinearStructs<las::MICRO_BACKEND> * vecs);

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
  // template <typename I>
  // void freeeRVEBC(I bnd_bgn, I bnd_end, apf::Numbering * num);
  // TODO Move this to the implicit analysis since this only makes sense
  // when working with implicit solves
  void applyGuessSolution(FiberRVEAnalysis * ans,
                          const DeformationGradient & dfmGrd);
}  // namespace bio
#include "bioFiberRVEAnalysis_impl.h"
#endif
