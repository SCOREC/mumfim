#include "MultiscaleTissueAnalysis.h"
#include <Solvers.h>
#include <amsiCasters.h>
#include <amsiMultiscale.h>
#include <apfFunctions.h>
#include <apfNumbering.h>
#include <apfWrapper.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include "AnalysisIO.h"
#include "ModelTraits.h"
#include "MultiscaleConvergence.h"
#include "MultiscaleTissue.h"
#include "VolumeConvergence.h"
namespace mumfim
{
  void MultiscaleTissueIteration::iterate()
  {
    if (!PCU_Comm_Self())
      std::cout << "Multiscale Nonlinear Iteration : " << iteration()
                << std::endl;
    if (iteration() == 0) tssu->updateMicro();
    las->iter();
    fem_iter->iterate();
    tssu->iter();
    amsi::Iteration::iterate();
  }
  MultiscaleTissueAnalysis::MultiscaleTissueAnalysis(
      apf::Mesh * mesh,
      std::unique_ptr<mt::CategoryNode> analysis_case,
      MPI_Comm cm,
      const amsi::Analysis & amsi_analysis,
      const amsi::Multiscale & amsi_multiscale)
      : TissueAnalysis(mesh, std::move(analysis_case), cm, amsi_analysis)
      , cplng(getRelationID(amsi_multiscale.getMultiscaleManager(),
                            amsi_multiscale.getScaleManager(),
                            "macro",
                            "micro_fo"))
      , multiscale_(amsi_multiscale)
  {
    const auto * solution_strategy =
        mt::GetPrimaryCategoryByType(analysis_case.get(), "solution strategy");
    tssu = new MultiscaleTissue(mesh, *analysis_case, cm, multiscale_);
    addVolumeTracking(mesh,solution_strategy);
    // compute the multiscale tissue iteration after the volumes have been
    // computed
    itr_stps.push_back(new MultiscaleTissueIteration(
        static_cast<MultiscaleTissue *>(tssu), las));
    // checkpoint after performing an iteration (this way numbering lines up
    // properly)
    itr_stps.push_back(new TissueCheckpointIteration(this));
    itr = new amsi::MultiIteration(itr_stps.begin(), itr_stps.end());
    buildLASConvergenceOperators(solution_strategy, itr, las,
                                 std::back_inserter(cvg_stps));
    buildVolConvergenceOperators(solution_strategy, itr, tssu->getUField(),
                                 trkd_vols, std::back_inserter(cvg_stps));
    cvg = new MultiscaleConvergence(cvg_stps.begin(), cvg_stps.end(), cplng);
    static_cast<MultiscaleTissue *>(tssu)->initMicro();
  }
  void MultiscaleTissueAnalysis::finalizeStep()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->scaleBroadcast(cplng, &completed);
  }
}  // namespace mumfim
