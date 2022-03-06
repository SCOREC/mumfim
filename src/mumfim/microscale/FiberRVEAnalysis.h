#ifndef MUMFIM_FIBER_RVE_ANALYSIS_H_
#define MUMFIM_FIBER_RVE_ANALYSIS_H_
#include <amsiNonlinearAnalysis.h>
#include <lasSparskit.h>
#include <lasSparskitExterns.h>
#include <string>
#include "FiberNetwork.h"
#include "MassIntegrator.h"
#include "MicroFOParams.h"
#include "RVE.h"
#include "TrussIntegrator.h"
#include "MultiscaleMicroFOParams.h"
#include "RVEAnalysis.h"
namespace mumfim
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
    LinearStructs(FiberNetwork* fiber_network,
                  int ndofs,
                  double solver_tol,
                  std::shared_ptr<void> bfrs = nullptr);
    las::Mat * getK() const {if (k == nullptr){std::cerr<<"K BROKE LS"<<std::endl; std::abort();} return k; }
    las::Vec * getU() const { return u; }
    las::Vec * getF() const { return f; }
    // swaps the vector pointers if zero is true the v2 vector
    // will be zero'd
    las::Solve * getSlv() const { return slv.get(); }
    ~LinearStructs();
    // give direct access to buffers to avoid too much function call overhead
    // friend class FiberRVEAnalysis;
    // friend class FiberRVEAnalysisSImplicit;
    protected:
    las::Mat * k;
    las::Vec * u;
    las::Vec * f;
    std::unique_ptr<las::Solve> slv;
    //las::Solve * slv;
    private:
    // sparsity pattern of the network
    // Eventually this may be promoted to a shared pointer
    // because sparsity patterns are the same for all networks of the same type
    using sparsity_ptr_type = std::unique_ptr<las::Sparsity, decltype(&las::destroySparsity<las::MICRO_BACKEND>)>;
    sparsity_ptr_type mSparsity;
    std::shared_ptr<void> mBuffers;
  };
  class FiberRVEAnalysis : public RVEAnalysis
  {
    protected:
    std::unique_ptr<FiberNetwork> mFiberNetwork;
    std::unique_ptr<MicroSolutionStrategy> mSolutionStrategy;
    std::unique_ptr<RVE> rve;
    virtual void computeCauchyStress(double sigma[6]);

    public:
    std::vector<apf::MeshEntity *> bnd_nds[RVE::side::all + 1];
    double solver_eps;
    double prev_itr_factor;
    int max_cut_attempt;
    int attempt_cut_factor;
    int max_itrs;
    amsi::DetectOscillationType detect_osc_type;
    FiberRVEAnalysis(const FiberRVEAnalysis & an);
    FiberRVEAnalysis(std::unique_ptr<FiberNetwork> fn,
                     std::unique_ptr<MicroSolutionStrategy> ss,
                     std::shared_ptr<void> workspace);
    virtual ~FiberRVEAnalysis();
    FiberRVEAnalysis& operator=(FiberRVEAnalysis&& other);
    FiberNetwork * getFn() const { return mFiberNetwork.get(); }
    double getCurrentVolume() const { return rve->measureDu(); }
    // TODO make this protected, makes the api more clear
    RVE * getRVE() const { return rve.get(); }
    virtual bool run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords=true) override = 0;
    virtual SolverType getAnalysisType() = 0;

    std::shared_ptr<LinearStructs<las::MICRO_BACKEND>> vecs;
    //LinearStructs<las::MICRO_BACKEND> * vecs;
    // get the global stiffness matrix
    las::Mat * getK() const { return vecs->getK(); }
    // get the global displacement vector
    las::Vec * getU() const { return vecs->getU(); }
    // get the global force vector
    las::Vec * getF() const { return vecs->getF(); }
    // solve is needed for computing derivs
    las::Solve * getSlv() const { return vecs->getSlv(); }
  };
  std::unique_ptr<FiberRVEAnalysis> createFiberRVEAnalysis(
      std::unique_ptr<FiberNetwork> fiber_network,
      std::unique_ptr<MicroSolutionStrategy> solution_strategy,
      std::shared_ptr<void> workspace = nullptr);
  std::unique_ptr<LinearStructs<las::MICRO_BACKEND>>
    createLinearStructs(
      FiberNetwork * fiber_network,
      int ndofs,
      double solver_tol,
      shared_ptr<void> bfrs = nullptr);
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
}  // namespace mumfim
#include "FiberRVEAnalysis_impl.h"
#endif
