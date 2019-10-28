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
  /*
   * Analysis type enum describes how the solver should function, e.g.
   * defaul explicit or implicit types. In the future this could select
   * between different subtypes as well e.g. ImplicitNewton and
   * ImplicitContinuation
   */
  enum class FiberRVEAnalysisType
  {
    StaticImplicit,      ///< Newton implicit solver
    Explicit  ///< Central difference explicit solver
  };
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
                  void * bfrs = NULL,
                  las::Sparsity * massCsr = NULL);
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
    FiberNetwork * fn;
    RVE * rve;
    virtual void computeCauchyStress(double sigma[6]);

    public:
    std::vector<apf::MeshEntity *> bnd_nds[RVE::side::all + 1];
    apf::Integrator * es;
    LinearStructs<las::MICRO_BACKEND> * vecs;
    double solver_eps;
    double prev_itr_factor;
    int max_cut_attempt;
    int attempt_cut_factor;
    int max_itrs;
    amsi::DetectOscillationType detect_osc_type;
    FiberRVEAnalysis(const FiberRVEAnalysis & an);
    FiberRVEAnalysis(FiberNetwork * fn,
                     LinearStructs<las::MICRO_BACKEND> * vecs,
                     const MicroSolutionStrategy & ss);
    virtual ~FiberRVEAnalysis();
    // get the global stiffness matrix
    las::Mat * getK() const { return vecs->k; }
    // get the global displacement vector
    las::Vec * getU() const { return vecs->u; }
    // get the global force vector
    las::Vec * getF() const { return vecs->f; }
    // solve is needed for computing derivs
    las::Solve * getSlv() const { return vecs->slv; }
    FiberNetwork * getFn() const { return fn; }
    RVE * getRVE() const { return rve; }
    virtual bool run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords=true) override = 0;
    virtual FiberRVEAnalysisType getAnalysisType() = 0;
    // these functions need to be implemented by any classes
    // that use alternative storage to the las vec, or don't explicitly
    // compute the stiffness matrix.
    virtual void computeStiffnessMatrix()
    {
      auto ops = las::getLASOps<las::MICRO_BACKEND>();
      ops->zero(getK());
      ops->zero(getF());
      es->process(getFn()->getNetworkMesh(), 1);
      // finalize the vectors so we can set boundary condition
      // values
      las::finalizeMatrix<las::MICRO_BACKEND>(getK());
      las::finalizeVector<las::MICRO_BACKEND>(getF());
    }
  };
  class FiberRVEAnalysisSImplicit : public FiberRVEAnalysis
  {
    protected:
    virtual void computeCauchyStress(double sigma[6]) final;

    public:
    // constructors
    explicit FiberRVEAnalysisSImplicit(const FiberRVEAnalysisSImplicit & an);
    FiberRVEAnalysisSImplicit(FiberNetwork * fn,
                              LinearStructs<las::MICRO_BACKEND> * vecs,
                              const MicroSolutionStrategy & ss);
    virtual bool run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords=true) final;
    virtual FiberRVEAnalysisType getAnalysisType()
    {
      return FiberRVEAnalysisType::StaticImplicit;
    }
  };
  FiberRVEAnalysis * createFiberRVEAnalysis(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      micro_fo_solver & slvr,
      micro_fo_int_solver & slvr_int);
  FiberRVEAnalysis * createFiberRVEAnalysis(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      const MicroSolutionStrategy & ss,
      FiberRVEAnalysisType type = FiberRVEAnalysisType::StaticImplicit);
  FiberRVEAnalysis * initFiberRVEAnalysisFromMultiscale(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      micro_fo_header & hdr,
      micro_fo_params & prm,
      micro_fo_init_data & ini,
      micro_fo_solver & slvr,
      micro_fo_int_solver & slvr_int);
  LinearStructs<las::MICRO_BACKEND> * createLinearStructs(
      int ndofs,
      double solver_tol,
      las::Sparsity * csr,
      void * bfrs = NULL,
      las::Sparsity * massSprs = NULL);
  void destroyFiberRVEAnalysisLinearStructs(
      LinearStructs<las::MICRO_BACKEND> * vecs);
  /*
   * perform a deep copy of the fiber rve analysis
   */
  // FiberRVEAnalysis * copyAnalysis(FiberRVEAnalysis * an);

  class FiberRVEIterationSImplicit : public amsi::Iteration
  {
    protected:
    FiberRVEAnalysisSImplicit * an;

    public:
    FiberRVEIterationSImplicit(FiberRVEAnalysisSImplicit * a);
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
  // template <typename I>
  // void freeeRVEBC(I bnd_bgn, I bnd_end, apf::Numbering * num);
  // TODO Move this to the implicit analysis since this only makes sense
  // when working with implicit solves
  void applyGuessSolution(FiberRVEAnalysis * ans,
                          const DeformationGradient & dfmGrd);
}  // namespace bio
#include "bioFiberRVEAnalysis_impl.h"
#endif
