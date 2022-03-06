#include "bioTissueAnalysis.h"
#include <apf.h>
#include <iostream>
#include "bioAnalysisIO.h"
#include "bioModelTraits.h"
namespace bio
{
  TissueAnalysis::TissueAnalysis(apf::Mesh * mesh,
                                 std::unique_ptr<const mt::CategoryNode> cs,
                                 MPI_Comm c)
      : cm(c)
      , mesh(mesh)
      , analysis_case(std::move(cs))
      , t(0.0)
      , dt(0.0)
      , stp(0)
      , mx_stp(1)
      , tssu(NULL)
      , itr()
      , itr_stps()
      , cvg()
      , cvg_stps()
      , trkd_vols()
      , las(new amsi::PetscLAS(0, 0))
      , completed(false)
      , state_fn()
      , constraint_fn(amsi::fs->getResultsDir() + "/constraints.log")
      , frcs_fn(amsi::fs->getResultsDir() + "/loads.log")
      , nrms_fn(amsi::fs->getResultsDir() + "/norms.log")
      , dsps_fn(amsi::fs->getResultsDir() + "/disps.log")
      , vols_fn(amsi::fs->getResultsDir() + "/vols.log")
      , state()
      , constraints()
      , frcs()
      , nrms()
      , dsps()
      , vols()
      , frc_itms()
      , dsp_itms()
      , vol_itms()
  {
  }
  TissueAnalysis::~TissueAnalysis()
  {
    // since we know all of the iteration steps are allocated on the heap delete
    // them
    for (auto itr_stp = itr_stps.begin(); itr_stp != itr_stps.end(); ++itr_stp)
    {
      delete (*itr_stp);
      (*itr_stp) = NULL;
    }
    delete itr;
    // since we know all of the convegence steps are allocated on the heap
    // delete them
    for (auto cvg_stp = cvg_stps.begin(); cvg_stp != cvg_stps.end(); ++cvg_stp)
    {
      delete (*cvg_stp);
      (*cvg_stp) = NULL;
    }
    delete cvg;
    delete tssu;
    delete las;
  }
  void TissueAnalysis::init()
  {
    // util data
    int rnk = -1;
    MPI_Comm_rank(cm, &rnk);
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
    tssu = new NonlinearTissue(mesh, *analysis_case, cm);
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
    // We want to do the tissue iteration after we compute the volumes
    itr_stps.push_back(new TissueIteration(tssu, las));
    itr_stps.push_back(new TissueCheckpointIteration(this));
    itr = new amsi::MultiIteration(itr_stps.begin(), itr_stps.end());
    buildLASConvergenceOperators(solution_strategy,itr,las,std::back_inserter(cvg_stps));
    buildVolConvergenceOperators(solution_strategy,itr,tssu->getUField(),trkd_vols,std::back_inserter(cvg_stps));
    cvg = new amsi::MultiConvergence(cvg_stps.begin(), cvg_stps.end());
    // output params
#ifdef LOGRUN
    std::stringstream cnvrt;
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
    amsi::log(state) << "STEP, ITER,   T, DESC" << std::endl
                     << "   0,    0, 0.0, init" << std::endl;
#endif
  }
  void TissueAnalysis::run()
  {
    tssu->preRun();
    tssu->recoverSecondaryVariables(stp);
    checkpoint();
    // write the initial state of everything
    t += dt;
    tssu->setSimulationTime(t);
    logVolumes(vol_itms.begin(), vol_itms.end(), vols, stp, tssu->getUField());
    tssu->computeInitGuess(las);
    completed = false;
    while (!completed)
    {
#ifdef LOGRUN
      amsi::log(state) << stp << ", " << itr->iteration() << ", " << MPI_Wtime()
                       << ", "
                       << "start_step" << std::endl;
#endif
      if (!PCU_Comm_Self()) std::cout << "Load step = " << stp << std::endl;
      if (amsi::numericalSolve(itr, cvg))
      {
        // checkpoint the initial state
        // note this is not actually the initial state
        // since we have already applied our guess solution
        // if(stp == 0 && itr->iteration() == 0)
#ifdef LOGRUN
        amsi::log(state) << stp << ", " << itr->iteration() << ", "
                         << MPI_Wtime() << ", "
                         << "end_solve" << std::endl;
#endif
        if (stp >= mx_stp - 1)
        {
          completed = true;
          std::cout << "Final load step converged. Case complete." << std::endl;
        }
        else
        {
          for (auto vol = trkd_vols.begin(); vol != trkd_vols.end(); ++vol)
            vol->second->step();
          las->step();
          tssu->step();
        }
#ifdef LOGRUN
        amsi::log(state) << stp << ", " << itr->iteration() << ", "
                         << MPI_Wtime() << ", "
                         << "end_step" << std::endl;
#endif
      }
      else
      {
        completed = true;
        std::cerr << "ERROR: Step " << stp << " failed to converge!"
                  << std::endl;
        finalizeStep();
      }
      logDisps(dsp_itms.begin(), dsp_itms.end(), dsps, stp, tssu->getUField());
      logForces(frc_itms.begin(), frc_itms.end(), frcs, stp, tssu);
      logVolumes(vol_itms.begin(), vol_itms.end(), vols, stp,
                 tssu->getUField());
      std::cout << "checkpointing (macro)" << std::endl;
      std::cout << "Rewriting at end of load step to include orientation data"
                << std::endl;
      tssu->recoverSecondaryVariables(stp);
      checkpoint();
      // reset the iteration from the numerical solve after checkpointing which
      // records iteration information
      itr->reset();
      stp++;
      t += dt;
      tssu->setSimulationTime(t);
      std::cout << "Finalizing step (macro)" << std::endl;
      finalizeStep();
    }
    deinit();
  }
  void TissueAnalysis::finalizeStep(){};
  void TissueAnalysis::checkpoint()
  {
#ifdef LOGRUN
    int rnk = -1;
    MPI_Comm_rank(cm, &rnk);
    if (rnk == 0)
    {
      std::ofstream frcs_fs(frcs_fn.c_str(), std::ios::out | std::ios::app);
      std::ofstream dsps_fs(dsps_fn.c_str(), std::ios::out | std::ios::app);
      std::ofstream vols_fs(vols_fn.c_str(), std::ios::out | std::ios::app);
      std::ofstream nrms_fs(nrms_fn.c_str(), std::ios::out | std::ios::app);
      std::ofstream cnst_fs(constraint_fn.c_str(),
                            std::ios::out | std::ios::app);
      amsi::flushToStream(frcs, frcs_fs);
      amsi::flushToStream(dsps, dsps_fs);
      amsi::flushToStream(vols, vols_fs);
      amsi::flushToStream(nrms, nrms_fs);
      amsi::flushToStream(constraints, cnst_fs);
    }
    std::ofstream st_fs(state_fn.c_str(), std::ios::out | std::ios::app);
    amsi::flushToStream(state, st_fs);
#endif
    // write mesh to file
    std::string pvd("/out.pvd");
    std::ofstream pvdf;
    int iteration = itr->iteration() - 1;
    std::cout << "ITERATION: " << iteration << std::endl;
    std::stringstream cnvrt;
    cnvrt << "msh_stp_" << stp << "_iter_";
    // amsi::writePvdFile(pvd, cnvrt.str(), iteration-1);
    cnvrt << iteration;
    apf::writeVtkFiles(
        std::string(amsi::fs->getResultsDir() + "/" + cnvrt.str()).c_str(),
        tssu->getMesh());
  }
  void TissueAnalysis::revert() {}
  void TissueAnalysis::deinit()
  {
#ifdef LOGRUN
    int rnk = -1;
    MPI_Comm_rank(cm, &rnk);
    if (rnk == 0)
    {
      amsi::deleteLog(vols);
      amsi::deleteLog(dsps);
      amsi::deleteLog(nrms);
      amsi::deleteLog(frcs);
      amsi::deleteLog(constraints);
    }
    amsi::deleteLog(state);
#endif
  }
}  // namespace bio
