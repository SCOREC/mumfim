#ifndef BIO_TISSUE_ANALYSIS_H_
#define BIO_TISSUE_ANALYSIS_H_
#include <amsiNonlinearAnalysis.h>
#include "bioNonlinearTissue.h"
#include "bioVolumeConvergence.h"
#include <model_traits/CategoryNode.h>
namespace bio
{
  /**
   * Extracts the volume convergence regions from a simmetrix case
   * @param solution_strategy associated solution strategy node
   * @param las the amsi las system used for convergence
   * @param out typically std::back_inserter to some sort of list type.
   */
  template <typename O>
  void buildLASConvergenceOperators(
      const mt::CategoryNode * solution_strategy,
      amsi::MultiIteration * it,
      amsi::LAS * las,
      O out);

  /**
   * Extracts the volume convergence regions from a simmetrix case
   * @param solution_strategy associated solution strategy node
   * @param it amsi iterator
   * @param u displacement
   * @param vl_tks list of volumes tracked for volume convergence
   * @param out typically std::back_inserter to some sort of list type.
   */
  template <typename I, typename O>
  void buildVolConvergenceOperators(
      const mt::CategoryNode * solution_strategy,
      amsi::MultiIteration * it,
      apf::Field * u,
      I vl_tks,
      O out);

  class TissueAnalysis
  {
    public:
    TissueAnalysis(apf::Mesh* mesh, std::unique_ptr<const mt::CategoryNode> analysis_case, MPI_Comm);
    ~TissueAnalysis();
    virtual void run();
    virtual void init();
    virtual void checkpoint();
    virtual void finalizeStep();
    virtual void revert();
    virtual void deinit();

    protected:
    // util
    MPI_Comm cm;
    std::unique_ptr<const mt::CategoryNode> analysis_case;
    apf::Mesh* mesh;
    // analysis
    double t;
    double dt;
    int stp;
    int mx_stp;
    NonlinearTissue * tssu;
    // this is actually a MultiIteration
    amsi::MultiIteration * itr;
    std::vector<amsi::Iteration *> itr_stps;
    // this is actually a MultiConvergence
    amsi::Convergence * cvg;
    std::vector<amsi::Convergence *> cvg_stps;
    // name of track volume model trait and ptr to volume calc
    std::map<std::string, VolCalc *> trkd_vols;
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
    std::vector<apf::ModelEntity *> frc_itms;
    std::vector<apf::ModelEntity *> dsp_itms;
    std::vector<apf::ModelEntity *> vol_itms;
  };
  class TissueIteration : public amsi::Iteration
  {
    protected:
    NonlinearTissue * tssu;
    amsi::LAS * las;

    public:
    TissueIteration(NonlinearTissue * t, amsi::LAS * l) : tssu(t), las(l) {}
    virtual void iterate()
    {
      LinearSolver(tssu, las);
      tssu->iter();
      las->iter();
      amsi::Iteration::iterate();
    }
  };
  class TissueCheckpointIteration : public amsi::Iteration
  {
    protected:
    TissueAnalysis * tssu;

    public:
    TissueCheckpointIteration(TissueAnalysis * t) : tssu(t) {}
    virtual void iterate()
    {
      std::cout << "Checkpointing iteration: " << this->iteration()
                << std::endl;
      tssu->checkpoint();
      amsi::Iteration::iterate();
    }
  };
}  // namespace bio
#include "bioTissueAnalysis_impl.h"
#endif
