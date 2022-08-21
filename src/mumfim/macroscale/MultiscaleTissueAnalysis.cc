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
#include "mumfim/exceptions.h"
namespace mumfim
{
  void MultiscaleTissueIteration::iterate()
  {
    std::cerr<<"We shouldn't be calling the multiscale tissue iteration if code is working as expected\n";
    std::abort();
    if (!PCU_Comm_Self())
      std::cout << "Multiscale Nonlinear Iteration : " << iteration()
                << std::endl;
    // Note!: We shouldn't need to call this here since preRun calls update micro!!
    if (iteration() == 0) tssu->updateMicro();
    // multiscale calls las->iter *before* fem_iter iterate
    // single scale calls it afet! does this matter?
    // This should only impact the first iteration. Where multiscale is probably
    // wrong. Thaat's because this is used to updat the vectors for convergence.
    las->iter();
    fem_iter->iterate();
    tssu->iter(); // ends up calling nonlineartissue iter which does volume stuff don't care about amsi::Iteration::iterate();
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
        mt::GetPrimaryCategoryByType((this->analysis_case.get()), "solution strategy");
    if(solution_strategy == nullptr) {
      throw mumfim_error("solution strategy is required for multiscale analysis");
    }
    tssu = new MultiscaleTissue(mesh, *(this->analysis_case), cm, multiscale_);
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

  void MultiscaleTissueAnalysis::finalizeIteration(bool lastIteration) {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->scaleBroadcast(cplng, &lastIteration);
  }
}  // namespace mumfim
