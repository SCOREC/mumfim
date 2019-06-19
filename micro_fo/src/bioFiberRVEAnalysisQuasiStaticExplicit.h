#ifndef BIO_FIBER_RVE_ANALYSIS_QUASI_STATIC_EXPLICIT
#define BIO_FIBER_RVE_ANALYSIS_QUASI_STATIC_EXPLICIT
#include "bioFiberRVEAnalysis.h"
namespace bio
{
  apf::Integrator * createExplicitMicroElementalSystem(FiberNetwork * fn,
                                                       las::Mat * k,
                                                       las::Vec * f,
                                                       las::Vec * f_damp,
                                                       las::Vec * v);
  class FiberRVEAnalysisQSExplicit : public FiberRVEAnalysis
  {
    protected:
    double internal_energy;
    double kinetic_energy;
    double external_energy;  // applied work
    double damping_energy;
    double time;  // t(n)
    double stiffness_damping_factor;
    double mass_damping_factor;
    las::Vec * m_inv;
    public:
    // explicit FiberRVEAnalysisQSExplicit(const FiberRVEAnalysisSImplicit &
    // an);
    FiberRVEAnalysisQSExplicit(FiberNetwork * fn,
                               LinearStructs<las::MICRO_BACKEND> * vecs,
                               const MicroSolutionStrategy & ss);
    virtual ~FiberRVEAnalysisQSExplicit();
    // get the global mass matrix
    las::Mat * getM() const { return vecs->m; }
    las::Vec * getMInv() const { return m_inv; }
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
    // for testing purposes
    protected:
    apf::Field * f_int_field;
    apf::Field * f_ext_field;
    apf::Field * f_damp_field;
    apf::Field * f_inertial_field;
    las::Vec * f_inertial;

    public:
    apf::Field * getFIntField() { return f_int_field; }
    apf::Field * getFExtField() { return f_ext_field; }
    apf::Field * getFDampField() { return f_damp_field; }
    apf::Field * getFInertialField() { return f_inertial_field; }
    las::Vec * getFInertial() { return f_inertial; }
    void updateFInt()
    {
      vecs->swapVec(vecs->f_int, vecs->prev_f_int);
      static_cast<ExplicitTrussIntegrator *>(es)->updateF(getFInt());
    }
    void updateFExt() { vecs->swapVec(vecs->f_ext, vecs->prev_f_ext); }
    void updateFDamp()
    {
      vecs->swapVec(vecs->f_damp, vecs->prev_f_damp);
      static_cast<ExplicitTrussIntegrator *>(es)->updateFDamp(getFDamp());
    }
    double getTotalEnergy() const
    {
      return internal_energy + kinetic_energy + external_energy;
    }
    virtual FiberRVEAnalysisType getAnalysisType()
    {
      return FiberRVEAnalysisType::QuasiStaticExplicit;
    }
  };
  class ExplicitOutputWriter;
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
    unsigned int history_print_steps;
    unsigned int field_print_steps;
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
                                unsigned int historyPrintSteps,
                                unsigned int fieldPrintSteps,
                                ExplicitOutputWriter * outputWriter);
    virtual void iterate();
    bool getCompleted() { return last_iter; }
  };
  class ExplicitOutputWriter
  {
    private:
    unsigned int outputFrame;
    std::vector<amsi::PvdData> pvdData;
    std::string folder;
    std::string pvdName;
    FiberRVEAnalysisQSExplicit * an;
    las::LasOps<las::MICRO_BACKEND> * ops;

    public:
    ExplicitOutputWriter(std::string folder,
                         std::string pvdName,
                         FiberRVEAnalysisQSExplicit * an);
    void writeHistoryData(int iteration);
    void writeFieldData(int iteration);
  };
  template <typename I>
  void applyBndNdsToVec(apf::Numbering * num,
                        apf::Field * field,
                        las::Vec * vec,
                        I bnd_bgn,
                        I bnd_end)
  {
    /*
    apf::Mesh * mesh = apf::getMesh(field);
    apf::FieldShape * shp = apf::getShape(field);
    int ncomp = apf::countComponents(field);
    // double * cmps = new double[ncomp];
    double cmps[ncomp];
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    int dof;
    for (I ent = bnd_bgn; ent != bnd_end; ++ent)
    {
      int nnodes = shp->countNodesOn(mesh->getType(*ent));
      for (int nd = 0; nd < nnodes; ++nd)
      {
        apf::getComponents(field, *ent, nd, cmps);
        for (int cmp = 0; cmp < ncomp; ++cmp)
        {
          if (apf::isNumbered(num, *ent, nd, cmp))
          {
            dof = apf::getNumber(num, *ent, nd, cmp);
            ops->set(vec, 1, &dof, &(cmps[cmp]));
          }
        }
      }
    }
    */
  }
  class ApplyDeformationGradientExplicit : public ApplyDeformationGradient
  {
    protected:
    apf::Matrix3x3 dVdX;
    apf::Matrix3x3 dAdX;
    apf::Field * vField;
    apf::Field * aField;
    int numComp;

    public:
    template <class IteratorType>
    ApplyDeformationGradientExplicit(IteratorType ent_bgn,
                                     IteratorType ent_end,
                                     apf::Matrix3x3 F,
                                     Amplitude * amp,
                                     double time,
                                     apf::Mesh * mesh,
                                     apf::Field * uField,
                                     apf::Field * duField,
                                     apf::Field * vField,
                                     apf::Field * aField)
        : ApplyDeformationGradient(ent_bgn, ent_end, F, mesh, duField, uField)
        , vField(vField)
        , aField(aField)
    {
      // only apply the portion of the deformation that we expect at time t
      FmI = FmI * (*amp)(time);
      dVdX = FmI * amp->derivative(time);
      dAdX = FmI * amp->secondDerivative(time);
    }
    virtual void atNode(int nd)
    {
      ApplyDeformationGradient::atNode(nd);
      // apply velocity and acceleration fields
      apf::Vector3 nd_v = dVdX * nd_xyz;
      apf::Vector3 nd_a = dAdX * nd_xyz;
      apf::setVector(vField, ent, nd, nd_v);
      apf::setVector(aField, ent, nd, nd_a);
    }
  };
}  // namespace bio
#endif
