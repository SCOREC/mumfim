#ifndef BIO_TISSUE_ANALYSIS_H_
#define BIO_TISSUE_ANALYSIS_H_
#include "bioNonlinearTissue.h"
#include "bioVolumeConvergence.h"
#include <amsiNonlinearAnalysis.h>
namespace bio
{
  template <typename O>
    void buildLASConvergenceOperators(pACase ss, amsi::Iteration * it, amsi::LAS * las, O out);
  template <typename I, typename O>
    void buildVolumeConvergenceOperators(pACase ss, amsi::Iteration * it, I vl_tks, apf::Field * u, O out);
  class TissueIteration : public amsi::ModularIteration
  {
  protected:
    NonlinearTissue * tssu;
    amsi::LAS * las;
  public:
    TissueIteration(NonlinearTissue * t, amsi::LAS * l)
      : amsi::ModularIteration()
      , tssu(t)
      , las(l)
    {}
    virtual void iterate()
    {
      LinearSolver(tssu,las);
      tssu->iter();
      las->iter();
      amsi::ModularIteration::iterate();
    }
  };
  class TissueAnalysis
  {
  public:
    TissueAnalysis(pGModel, pParMesh, pACase, MPI_Comm);
    ~TissueAnalysis();
    virtual void run();
    virtual void init();
    virtual void checkpoint();
    virtual void revert();
    virtual void deinit();
  protected:
    // util
    MPI_Comm cm;
    // simmetrix
    pGModel mdl;
    pParMesh msh;
    pACase cs;
    // analysis
    double t;
    double dt;
    int stp;
    int mx_stp;
    NonlinearTissue * tssu;
    amsi::Iteration * itr;
    amsi::Convergence * cvg;
    std::vector<amsi::Convergence*> cvg_stps;
    std::map<pANode,VolCalc*> trkd_vols;
    amsi::LAS * las;
    bool completed;
    // log filenames
    std::string state_fn;
    std::string constraint_fn;
    std::string frcs_fn;
    std::string nrms_fn;
    std::string dsps_fn;
    std::string vols_fn;
    // logs
    amsi::Log state;
    amsi::Log constraints;
    amsi::Log frcs;
    amsi::Log nrms;
    amsi::Log dsps;
    amsi::Log vols;
    // track model ents
    std::vector<apf::ModelEntity*> frc_itms;
    std::vector<apf::ModelEntity*> dsp_itms;
    std::vector<apf::ModelEntity*> vol_itms;
  };
}
#include "bioTissueAnalysis_impl.h"
#endif