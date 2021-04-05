#include "bioTissueMultiscaleAnalysis.h"
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
#include "bioAnalysisIO.h"
#include "bioModelTraits.h"
#include "bioMultiscaleConvergence.h"
#include "bioMultiscaleTissue.h"
#include "bioVolumeConvergence.h"
namespace bio
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
      MPI_Comm cm)
      : TissueAnalysis(mesh, std::move(analysis_case), cm)
      , cplng(getRelationID(amsi::getMultiscaleManager(),
                            amsi::getScaleManager(),
                            "macro",
                            "micro_fo"))
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
    tssu = new MultiscaleTissue(mesh, *analysis_case, cm);
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
    // output params
#ifdef LOGRUN
    std::ostringstream cnvrt;
    cnvrt << rnk;
    state_fn =
        amsi::fs->getResultsDir() + "/tissue_state." + cnvrt.str() + ".log";
    amsi::getTrackedModelItems(cs, "output force",
                               std::back_inserter(frc_itms));
    amsi::getTrackedModelItems(cs, "output displacement",
                               std::back_inserter(dsp_itms));
    amsi::getTrackedModelItems(cs, "output volume",
                               std::back_inserter(vol_itms));
    // initialize logging
    state = amsi::activateLog("tissue_efficiency");
    if (rnk == 0)
    {
      constraints = amsi::activateLog("constraints");
      frcs = amsi::activateLog("loads");
      nrms = amsi::activateLog("norms");
      dsps = amsi::activateLog("displacement");
      vols = amsi::activateLog("volume");
      amsi::log(constraints) << "STEP, ITER, LAMBDA, BETA" << std::endl;
      amsi::log(frcs) << "STEP, ENT, I, J, K" << std::endl;
      amsi::log(nrms) << "STEP, ENT, NRM" << std::endl;
      amsi::log(dsps) << "STEP, ENT, X, Y, Z" << std::endl;
      amsi::log(vols) << "STEP, ENT, VOL" << std::endl;
    }
    amsi::log(state) << "STEP, ITER,   T, DESC" << std::endl;
#endif
  }
  void MultiscaleTissueAnalysis::finalizeStep()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->scaleBroadcast(cplng, &completed);
  }
  void MultiscaleTissueAnalysis::run() { TissueAnalysis::run(); }
}  // namespace bio
