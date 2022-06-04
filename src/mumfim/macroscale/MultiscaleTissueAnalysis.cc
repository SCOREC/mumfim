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
  }
  // TODO move to constructor
  void MultiscaleTissueAnalysis::init()
  {
    const auto * problem_definition =
        mt::GetPrimaryCategoryByType(analysis_case.get(), "problem definition");
    const auto * solution_strategy =
        mt::GetPrimaryCategoryByType(analysis_case.get(), "solution strategy");
    if (problem_definition == nullptr || solution_strategy == nullptr ||
        problem_definition->GetType() != "macro" ||
        solution_strategy->GetType() != "macro")
    {
      std::cerr << "Analysis case should have  \"problem definition\" and "
                   "\"solution strategy\" of the \"macro\" analysis type.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    // analysis params
    tssu = new MultiscaleTissue(mesh, *analysis_case, cm, multiscale_);
    const auto * timesteps_trait = mt::GetCategoryModelTraitByType<mt::IntMT>(
        solution_strategy, "num timesteps");
    if (timesteps_trait == nullptr)
    {
      std::cerr << "\"solution strategy\" must have \"num timesteps\" trait";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    mx_stp = (*timesteps_trait)();
    dt = (double)1.0 / (double)mx_stp;
    const auto * track_volume =
        mt::GetCategoryByType(solution_strategy, "track volume");
    if (track_volume != nullptr)
    {
      for (const auto & tracked_volume : track_volume->GetModelTraitNodes())
      {
        std::vector<apf::ModelEntity *> model_entities;
        GetModelTraitNodeGeometry(mesh, &tracked_volume, model_entities);
        trkd_vols[tracked_volume.GetName()] = new VolCalc(
            model_entities.begin(), model_entities.end(), tssu->getUField());
      }
    }
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
