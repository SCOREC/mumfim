#include "mumfim/macroscale/TissueAnalysis.h"
#include <apf.h>
#include <iostream>
#include "mumfim/macroscale/AnalysisIO.h"
#include "mumfim/macroscale/ModelTraits.h"
namespace mumfim
{
  TissueAnalysis::TissueAnalysis(apf::Mesh * mesh,
                                 std::unique_ptr<const mt::CategoryNode> cs,
                                 MPI_Comm c,
                                 const amsi::Analysis & amsi_analysis)
      : cm(c)
      , analysis_case(std::move(cs))
      , mesh(mesh)
      , t(0.0)
      , dt(0.0)
      , stp(0)
      , mx_stp(1)
      , tssu(nullptr)
      , itr()
      , itr_stps()
      , cvg()
      , cvg_stps()
      , trkd_vols()
      , las(new amsi::PetscLAS(0, 0))
      , completed(false)
      , state_fn()
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
    const auto * timesteps_trait = mt::GetCategoryModelTraitByType<mt::IntMT>(
        solution_strategy, "num timesteps");
    if (timesteps_trait == nullptr)
    {
      std::cerr << R"("solution strategy" must have "num timesteps" trait)";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    mx_stp = (*timesteps_trait)();
    dt = (double)1.0 / (double)mx_stp;
    // output params
#ifdef LOGRUN
    int rnk = -1;
    std::stringstream cnvrt;
    MPI_Comm_rank(cm, &rnk);
    cnvrt << rnk;
    state_fn =
        amsi::fs->getResultsDir() + "/tissue_state." + cnvrt.str() + ".log";
    state = amsi::activateLog("tissue_efficiency");
    amsi::log(state) << "STEP, ITER,   T, DESC\n"<< "   0,    0, 0.0, init\n";
#endif
  }
  void TissueAnalysis::addVolumeTracking(
      apf::Mesh * mesh,
      const mt::CategoryNode * solution_strategy)
  {
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
  }
  TissueAnalysis::~TissueAnalysis()
  {
    // since we know all of the iteration steps are allocated on the heap delete
    // them
    for (auto itr_stp = itr_stps.begin(); itr_stp != itr_stps.end(); ++itr_stp)
    {
      delete (*itr_stp);
      (*itr_stp) = nullptr;
    }
    delete itr;
    // since we know all of the convegence steps are allocated on the heap
    // delete them
    for (auto cvg_stp = cvg_stps.begin(); cvg_stp != cvg_stps.end(); ++cvg_stp)
    {
      delete (*cvg_stp);
      (*cvg_stp) = nullptr;
    }
    delete cvg;
    delete tssu;
    delete las;
#ifdef LOGRUN
    amsi::deleteLog(state);
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
    //logVolumes(vol_itms.begin(), vol_itms.end(), vols, stp, tssu->getUField());
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
      //logDisps(dsp_itms.begin(), dsp_itms.end(), dsps, stp, tssu->getUField());
      //logForces(frc_itms.begin(), frc_itms.end(), frcs, stp, tssu);
      //logVolumes(vol_itms.begin(), vol_itms.end(), vols, stp,
      //           tssu->getUField());
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
  }
  void TissueAnalysis::finalizeStep(){};
  void TissueAnalysis::checkpoint()
  {
#ifdef LOGRUN
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
}  // namespace mumfim
