#ifndef MUMFIM_TISSUE_ANALYSIS_H_
#define MUMFIM_TISSUE_ANALYSIS_H_
#include <amsiNonlinearAnalysis.h>
#include "NonlinearTissue.h"
#include "VolumeConvergence.h"
#include <model_traits/CategoryNode.h>
namespace mumfim
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
    TissueAnalysis(apf::Mesh * mesh,
                   std::unique_ptr<const mt::CategoryNode> cs,
                   MPI_Comm c,
                   const amsi::Analysis & amsi_analysis);
    virtual ~TissueAnalysis();
    virtual void run();
    virtual void checkpoint();
    virtual void finalizeStep();
    virtual void finalizeIteration(int);

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
    int iteration{0};
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
    // logs
    amsi::Log state;

    void addVolumeTracking(apf::Mesh * mesh,
                           const mt::CategoryNode * solution_strategy);
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
      // Note this is same as FEMLinearIteration
      LinearSolver(tssu, las);

      tssu->iter();
      // copies LHS/RHS into "old" state as needed for computing convergence
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
}  // namespace mumfim
#include "TissueAnalysis_impl.h"
#endif
