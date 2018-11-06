#ifndef BIO_FIBER_RVE_ANALYSIS_H_
#define BIO_FIBER_RVE_ANALYSIS_H_
#include <amsiNonlinearAnalysis.h>
#include <lasSparskit.h>
#include <string>
#include "bioFiberNetwork.h"
#include "bioMassIntegrator.h"
#include "bioMicroFOParams.h"
#include "bioMultiscaleRVE.h"
#include "bioRVE.h"
#include "bioTrussIntegrator.h"
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
    QuasiStaticExplicit  ///< Central difference explicit solver
  };
  /* \brief create the elemental system integrator
   */
  apf::Integrator * createImplicitMicroElementalSystem(
      FiberNetwork * fn,
      las::Mat * k,
      las::Vec * f);
  apf::Integrator * createExplicitMicroElementalSystem(
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
    las::Mat * getM() const { return m; }
    las::Mat * getC() const { return c; }
    las::Vec * getU() const { return u; }
    las::Vec * getF() const { return f; }
    // swaps the vector pointers if zero is true the v2 vector
    // will be zero'd
    void swapVec(las::Vec *& v1, las::Vec *& v2, bool zero=false);
    las::Solve * getSlv() const { return slv; }
    ~LinearStructs();
    // give direct access to buffers to avoid too much function call overhead
    friend class FiberRVEAnalysis;
    friend class FiberRVEAnalysisSImplicit;
    friend class FiberRVEAnalysisQSExplicit;
    protected:
    las::Mat * k;
    las::Mat * m;
    las::Mat * c;
    las::Vec * u;
    las::Vec * f;
    las::Vec * delta_u;
    //las::Vec * prev_f;
    las::Vec * f_int;
    las::Vec * prev_f_int;
    las::Vec * f_ext;
    las::Vec * prev_f_ext;
    las::Vec * f_damp;
    las::Vec * prev_f_damp;
    las::Vec * v;
    las::Vec * prev_v;
    las::Vec * a;
    las::Solve * slv;
  };
  class FiberRVEAnalysis
  {
    protected:
    FiberNetwork * fn;

    public:
    // move this out of here! requires coupling data structures as currently
    // written hmmm...
    MultiscaleRVE * multi;
    RVE * rve;
    std::vector<apf::MeshEntity *> bnd_nds[RVE::side::all + 1];
    apf::Integrator * es;
    apf::DynamicMatrix dx_fn_dx_rve;
    LinearStructs<las::MICRO_BACKEND> * vecs;
    bool dx_fn_dx_rve_set;
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
    virtual bool run(const DeformationGradient & dfmGrd) = 0;
    virtual FiberRVEAnalysisType getAnalysisType() = 0;
  };
  class FiberRVEAnalysisSImplicit : public FiberRVEAnalysis
  {
    protected:
    public:
    // constructors
    explicit FiberRVEAnalysisSImplicit(const FiberRVEAnalysisSImplicit & an);
    FiberRVEAnalysisSImplicit(FiberNetwork * fn,
                              LinearStructs<las::MICRO_BACKEND> * vecs,
                              const MicroSolutionStrategy & ss);
    virtual bool run(const DeformationGradient & dfmGrd);
    virtual FiberRVEAnalysisType getAnalysisType()
    {
      return FiberRVEAnalysisType::StaticImplicit;
    }
  };
  class FiberRVEAnalysisQSExplicit : public FiberRVEAnalysis
  {
    protected:
    double internal_energy;
    double kinetic_energy;
    double external_energy;  // applied work
    double damping_energy;
    double time;             // t(n)
    public:
    // explicit FiberRVEAnalysisQSExplicit(const FiberRVEAnalysisSImplicit &
    // an);
    FiberRVEAnalysisQSExplicit(FiberNetwork * fn,
                               LinearStructs<las::MICRO_BACKEND> * vecs,
                               const MicroSolutionStrategy & ss);
    // get the global mass matrix
    las::Mat * getM() const { return vecs->m; }
    // get the velocity damping matrix
    las::Mat * getC() const { return vecs->c; }
    las::Vec * getA() const { return vecs->a; }
    las::Vec * getV() const { return vecs->v; }
    las::Vec * getFExt() const { return vecs->f_ext; }
    las::Vec * getPrevFExt() const { return vecs->prev_f_ext; }
    las::Vec * getFInt() const { return vecs->f_int; }
    las::Vec * getPrevFInt() const { return vecs->prev_f_int; }
    las::Vec * getFDamp() { return vecs->f_damp; }
    las::Vec * getPrevFDamp() { return vecs->prev_f_damp; }
    las::Vec * getDeltaU() { return vecs->delta_u; }
    void setC(las::Mat * c) { vecs->c = c; }
    virtual bool run(const DeformationGradient & dfmGrd);
    double getTime() const { return time; }
    void setTime(double t) { time = t; }
    double getInternalEnergy() const { return internal_energy; }
    double getKineticEnergy() const { return kinetic_energy; }
    double getExternalEnergy() const { return external_energy; }
    double getDampingEnergy() const { return damping_energy; }
    void setInternalEnergy(double e) { internal_energy = e; }
    void setKineticEnergy(double e) { kinetic_energy = e; }
    void setExternalEnergy(double e) { external_energy = e; }
    void setDampingEnergy(double e) { damping_energy = e; }
    void updateFInt()
    {
      vecs->swapVec(vecs->f_int, vecs->prev_f_int);
      static_cast<ExplicitTrussIntegrator *>(es)->updateF(getFInt());
    }
    void updateFExt() {
                        vecs->swapVec(vecs->f_ext, vecs->prev_f_ext);
                         }
    void updateFDamp() { vecs->swapVec(vecs->f_damp, vecs->prev_f_damp); }

    double getTotalEnergy() const
    {
      return internal_energy + kinetic_energy + external_energy;
    }
    virtual FiberRVEAnalysisType getAnalysisType()
    {
      return FiberRVEAnalysisType::QuasiStaticExplicit;
    }
  };
  FiberRVEAnalysis * createFiberRVEAnalysis(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      micro_fo_solver & slvr,
      micro_fo_int_solver & slvr_int,
      FiberRVEAnalysisType type = FiberRVEAnalysisType::StaticImplicit);
  FiberRVEAnalysis * createFiberRVEAnalysis(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      const MicroSolutionStrategy & ss,
      FiberRVEAnalysisType type = FiberRVEAnalysisType::StaticImplicit);
  FiberRVEAnalysis * initFromMultiscale(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      micro_fo_header & hdr,
      micro_fo_params & prm,
      micro_fo_init_data & ini,
      micro_fo_solver & slvr,
      micro_fo_int_solver & slvr_int);
  void destroyAnalysis(FiberRVEAnalysis *);
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
  class ExplicitOutputWriter;

  class FiberRVEIterationSImplicit : public amsi::Iteration
  {
    protected:
    FiberRVEAnalysisSImplicit * an;

    public:
    FiberRVEIterationSImplicit(FiberRVEAnalysisSImplicit * a);
    void iterate();
  };
  struct Amplitude
  {
    virtual double operator()(double time) = 0;
    virtual double derivative(double time) = 0;
    virtual double secondDerivative(double time) = 0;
    virtual double integral(double time) = 0;
    virtual ~Amplitude() {}
  };
  class FiberRVEIterationQSExplicit : public amsi::Iteration
  {
    private:
    las::Vec * tmp;
    protected:
    FiberRVEAnalysisQSExplicit * an;
    DeformationGradient appliedDefm;
    Amplitude * amplitude;
    bool last_iter;
    // print data every print_steps steps
    unsigned int print_steps;
    // the length of time the loading occurs for
    double load_time;
    // should be between 0.8 and 0.95 (Belytscheko)
    double time_step_factor;
    // Rayleigh damping factors
    double mass_damping_factor;
    double stiffness_damping_factor;
    ExplicitOutputWriter * outputWriter;
    public:
    FiberRVEIterationQSExplicit(FiberRVEAnalysisQSExplicit * a,
                                DeformationGradient appliedDefm,
                                Amplitude * amplitude,
                                double loadTime,
                                double timeStepFactor,
                                double massDampingFactor,
                                double stiffnessDampingFactor,
                                unsigned int printSteps,
                                ExplicitOutputWriter * outputWriter);
    void iterate();
    bool getCompleted() { return last_iter; }
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
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma);
  void applyGuessSolution(FiberRVEAnalysis * ans,
                          const DeformationGradient & dfmGrd);
}  // namespace bio
#include "bioFiberRVEAnalysis_impl.h"
#endif
