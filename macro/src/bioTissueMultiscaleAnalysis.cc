#include "bioTissueMultiscaleAnalysis.h"
#include "bioMultiscaleTissue.h"
#include <Solvers.h>
#include <ConvenienceFunctions.h>
#include <amsiMultiscale.h>
#include <apfWrapper.h>
#include <simAttributes.h>
#include <apfNumbering.h>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <fstream>
namespace bio
{
  TissueMultiScaleAnalysis::TissueMultiScaleAnalysis(pGModel imdl,pParMesh imsh,pACase cs,MPI_Comm cm)
    : rnk(-1)
    , rnkstr()
    , num_load_steps(1)
    , current_step(0)
    , t(0.0)
    , frc_itms()
    , dsp_itms()
    , vol_itms()
    , state(NULL)
    , cnstrnts(NULL)
    , norms(NULL)
    , disps(NULL)
    , loads(NULL)
    , vols(NULL)
    , state_file()
    , cnstrnts_file(amsi::fs->getResultsDir() + "/constraints.log")
    , norms_file(amsi::fs->getResultsDir() + "/norms.log")
    , disps_file(amsi::fs->getResultsDir() + "/disps.log")
    , loads_file(amsi::fs->getResultsDir() + "/loads.log")
    , vols_file(amsi::fs->getResultsDir() + "/vols.log")
    , model(imdl)
    , mesh(imsh)
    , part(PM_mesh(mesh,0))
    , tissue(NULL)
    , las(new amsi::PetscLAS(0,0))
  {
    MPI_Comm_rank(cm,&rnk);
    std::stringstream rnkstrm;
    rnkstrm << rnk;
    rnkstr = rnkstrm.str();
    state_file = amsi::fs->getResultsDir() + "/macro_state." + rnkstr + ".log";
    pACase pd = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
    tissue = new MultiscaleTissue(imdl,imsh,pd,cm);
    amsi::getTrackedModelItems(cs,"output force",std::back_inserter(frc_itms));
    amsi::getTrackedModelItems(cs,"output displacement",std::back_inserter(dsp_itms));
    amsi::getTrackedModelItems(cs,"output volume",std::back_inserter(vol_itms));
    pACase ss = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::SOLUTION_STRATEGY));
    num_load_steps = AttInfoInt_value((pAttInfoInt)AttNode_childByType((pANode)ss,"num timesteps"));
  }
  void TissueMultiScaleAnalysis::initLogs()
  {
    state = amsi::activateLog("macro_efficiency");
    amsi::log(state) << "TIMESTEP, ITERATION, TIME, STATUS, DESCRIPTION" << std::endl
                     << "0, 0, 0.0, ACTIVE, START" << std::endl;
    if(rnk == 0)
    {
      disps = amsi::activateLog("displacement");
      loads = amsi::activateLog("load");
      vols = amsi::activateLog("volume");
      norms = amsi::activateLog("norm_history");
      cnstrnts = amsi::activateLog("constraints");
      amsi::log(disps) << "LOADSTEP, ENT_TAG, X, Y, Z" << std::endl;
      amsi::log(loads) << "LOADSTEP, ENT_TAG, X, Y, Z" << std::endl;
      amsi::log(vols) << "LOADSTEP, ITERATION, ENT_TAG, INIT_VOL, VOL, VOL - PREV_VOL, VOL - INIT_VOL" << std::endl;
      amsi::log(norms) << "LOADSTEP, ITERATION, NORM" << std::endl;
      amsi::log(cnstrnts) << "LOADSTEP, ITERATION, LAMBDA, BETA" << std::endl;
    }
  }
  void TissueMultiScaleAnalysis::deleteLogs()
  {
    if(rnk == 0)
    {
      amsi::deleteLog(disps);
      amsi::deleteLog(loads);
      amsi::deleteLog(vols);
      amsi::deleteLog(norms);
      amsi::deleteLog(cnstrnts);
    }
    amsi::deleteLog(state);
  }
  int TissueMultiScaleAnalysis::run()
  {
    int result = 0;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    size_t cplng = cs->CommRelation_GetID("macro","micro_fo");
#   ifdef LOGRUN
    initLogs();
#   endif
    std::vector<std::pair< double, double > > norm_hist;
    // Run the analysis
    current_step = 0;
    updateTime();
    tissue->setSimulationTime(t);
    tissue->computeInitGuess(las);
    //apf::writeVtkFiles("init_guess",tissue->getMesh());
    tissue->initMicro();
    // Calculate Volume after initial guess
    tissue->updateVolumes();
    int complete = false;
    while(!complete)
    {
      bool converged = false;
      double eps = 1e-4;
      // incremental volume constraint
      // embedded_cell_resize:
      //double eps_v = 4.0;
      // accm. volume constraint V - V0
      double eps_v = 1e-3;
      unsigned iteration = 0;
      tissue->updateMicro();
#     ifdef LOGRUN
      if (rnk == 0)
        amsi::log(loads) << current_step << ", ";
#     endif
      /// Create convergence objects.
      VolumeConvergenceAccm_Incrmt dv_convergence(tissue,eps_v);
      LASResidualConvergence convergence(las,eps);
      while(!converged)
      {
#       ifdef LOGRUN
        amsi::log(state) << current_step << ", " << iteration << ", "
                         << amsi::getElapsedTime(state) << ", ACTIVE, BEGIN_ITER" << std::endl;
#       endif
        std::cout << "Current load step = " << current_step << "." << std::endl;
        std::cout << "Current nonlinear iteration = " << iteration << "." << std::endl;
        LinearSolver(tissue,las);
        las->iter();
        converged = convergence.converged();
        convergence.log(current_step, iteration, rnk);
        tissue->logCnstrntParams(current_step, iteration, rnk);
        // if we've converged on displacement, check the volume convergence and update the vols
        converged = converged ? dv_convergence.converged() : false;
        dv_convergence.log(current_step, iteration, rnk);
        cs->couplingBroadcast(cplng,&converged);
        tissue->iter();
        iteration++;
#       ifdef LOGRUN
        amsi::log(state) << current_step << ", " << iteration << ", "
                         << amsi::getElapsedTime(state) << ", ACTIVE, END_ITER" << std::endl;
#       endif
      }
      tissue->updatePrevVolumes();
#     ifdef LOGRUN
      if(rnk == 0)
      {
        std::fstream fs(state_file.c_str(), std::ios::out | std::ios::app);
        amsi::flushToStream(state,fs);
      }
      for(auto mdl_ent = dsp_itms.begin(); mdl_ent != dsp_itms.end(); ++mdl_ent)
      {
        //amsi::PrintField(tissue->getUField(),std::cout).run();
        std::list<pEntity> ents;
        int dm = amsi::modelItemTypeDim(GEN_type((pGEntity)*mdl_ent));
        amsi::getClassifiedDimEnts(part,(pGEntity)*mdl_ent,0,dm,std::back_inserter(ents));
        double dsp[3] = {};
        amsi::getAvgFieldValue(tissue->getUField(),ents.begin(),ents.end(),&dsp[0]);
        if(rnk == 0)
          amsi::log(disps) << current_step << ", " << GEN_tag((pGEntity)*mdl_ent) << ", " << dsp[0] << ", " << dsp[1] << ", " << dsp[2] << std::endl;
      }
      if(rnk == 0)
      {
        std::fstream fs(disps_file.c_str(), std::ios::out | std::ios::app);
        amsi::flushToStream(disps,fs);
      }
      for(auto mdl_ent = frc_itms.begin(); mdl_ent != frc_itms.end(); ++mdl_ent)
      {
        double frc[3] = {};
        tissue->getLoadOn((pGEntity)*mdl_ent,&frc[0]);
        if(rnk == 0)
          amsi::log(loads) << GEN_tag((pGEntity)*mdl_ent) << ", " << frc[0] << ", " << frc[1] << ", " << frc[2] << std::endl;
      }
      if(rnk == 0)
      {
        std::fstream lds_fs(loads_file.c_str(), std::ios::out | std::ios::app);
        std::fstream vls_fs(vols_file.c_str(), std::ios::out | std::ios::app);
        std::fstream nrms_fs(norms_file.c_str(), std::ios::out | std::ios::app);
        std::fstream cnst_fs(cnstrnts_file.c_str(), std::ios::out | std::ios::app);
        amsi::flushToStream(loads,lds_fs);
        amsi::flushToStream(vols,vls_fs);
        amsi::flushToStream(norms,nrms_fs);
        amsi::flushToStream(cnstrnts,cnst_fs);
      }
#     endif
      current_step++;
      // write mesh to file
      std::stringstream stpstrm;
      stpstrm << current_step;
      if ( (current_step) % 1 == 0 )
        apf::writeVtkFiles(std::string(amsi::fs->getResultsDir() + "/msh_stp_" + stpstrm.str()).c_str(),tissue->getMesh());
      if (current_step >= num_load_steps)
      {
        complete = true;
        std::cout << "Final load step converged! Simulation complete. Exiting..." << std::endl;
      }
      else
      {
        tissue->step();
        updateTime();
        tissue->setSimulationTime(t);
      }
      tissue->recoverSecondaryVariables(current_step);
      //      deleteLogs();
      // write mesh to file
      /*
        std::stringstream stpstrm;
        stpstrm << current_step;
        apf::writeVtkFiles(std::string(amsi::fs->getResultsDir() + "/msh_stp_" + stpstrm.str() + "_").c_str(),tissue->getMesh());
      */
      cs->couplingBroadcast(cplng,&complete);
    } // while(!complete)
    deleteLogs();
    return result;
  }
  void MultiscaleTissueIteration::iterate()
  {
    std::cout << "Nonlinear Iteration : " << iter << std::endl;
    LinearSolver(analysis->tissue_fea,analysis->linear_solver);
    // this really should be taken care of elsewhere...
    analysis->linear_solver->iter(); // accumulates the solution vector internally
    iter++;
  }
  bool MultiscaleTissueConvergence::converged()
  {
    double current_norm = 0.0;
    double accum_norm = 0.0;
    analysis->linear_solver->GetSolutionNorm(current_norm);
    analysis->linear_solver->GetAccumSolutionNorm(current_norm);
    int iteration = iter->getIteration();
    if(iteration >= 20 && iteration % 5 == 0)
      eps *= 10.0;
    bool converged = (current_norm / accum_norm) < eps;
    AMSI_DEBUG
      (
        std::cout << "Convergence condition: " << std::endl
        << "\t" << current_norm / accum_norm << " < " << eps << std::endl
        << "\t" << (converged ? "TRUE" : "FALSE") << std::endl;
        )
      amsi::ControlService * cs = amsi::ControlService::Instance();
    size_t cplng = cs->CommRelation_GetID("macro","micro_fo");
    cs->couplingBroadcast(cplng,&converged);
    return converged;
  }
} // end of namespace Biotissue
