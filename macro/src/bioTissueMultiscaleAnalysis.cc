#include "bioTissueMultiscaleAnalysis.h"
#include "bioAnalysisIO.h"
#include "bioMultiscaleTissue.h"
#include "bioVolumeConvergence.h"
#include <Solvers.h>
#include <ConvenienceFunctions.h>
#include <amsiCasters.h>
#include <amsiMultiscale.h>
#include <apfFunctions.h>
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
  void MultiscaleTissueIteration::iterate()
  {
    std::cout << "Multiscale Nonlinear Iteration : " << iteration() << std::endl;
    tssu->updateMicro();
    las->iter();
    LinearSolver(tssu,las);
    tssu->iter();
    amsi::Iteration::iterate();
  }
  TissueMultiScaleAnalysis::TissueMultiScaleAnalysis(pGModel imdl,pParMesh imsh,pACase cs,MPI_Comm cm)
    : rnk(-1)
    , num_load_steps(1)
    , current_step(0)
    , t(0.0)
    , frc_itms()
    , dsp_itms()
    , vol_itms()
    , state()
    , cnstrnts()
    , norms()
    , disps()
    , loads()
    , vols()
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
    state_file = amsi::fs->getResultsDir() + "/macro_state." + rnkstrm.str() + ".log";
    pACase pd = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
    tissue = new MultiscaleTissue(imdl,imsh,pd,cm);
    amsi::getTrackedModelItems(cs,"output force",std::back_inserter(frc_itms));
    std::vector<pModelItem> dp_itms;
    amsi::getTrackedModelItems(cs,"output displacement",std::back_inserter(dp_itms));
    std::transform(dp_itms.begin(),dp_itms.end(),std::back_inserter(dsp_itms),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    std::vector<pModelItem> vl_itms;
    amsi::getTrackedModelItems(cs,"output volume",std::back_inserter(vl_itms));
    std::transform(vl_itms.begin(),vl_itms.end(),std::back_inserter(vol_itms),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    pACase ss = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::SOLUTION_STRATEGY));
    num_load_steps = AttInfoInt_value((pAttInfoInt)AttNode_childByType((pANode)ss,"num timesteps"));
    itr = new MultiscaleTissueIteration(tissue,las);
    std::vector<pANode> cnvrg_nds;
    amsi::cutPaste<pANode>(AttNode_childrenByType((pANode)ss,"convergence operator"),std::back_inserter(cnvrg_nds));
    for(auto cnvrg_nd = cnvrg_nds.begin(); cnvrg_nd != cnvrg_nds.end(); ++cnvrg_nd)
    {
      char * tp = AttNode_imageClass(*cnvrg_nd);
      amsi::Convergence * cvg = NULL;
      if(std::string("volume convergence").compare(tp) == 0)
        cvg = buildBioConvergenceOperator(ss,*cnvrg_nd,itr,tissue->getUField());
      else
        cvg = amsi::buildSimConvergenceOperator(ss,*cnvrg_nd,itr,las);
      cnvrgs.push_back(cvg);
      Sim_deleteString(tp);
    }
    mlti_cvg = new amsi::MultiConvergence(cnvrgs.begin(),cnvrgs.end());
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
      amsi::log(vols) << "LOADSTEP, ENT_TAG, VOL" << std::endl;
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
  void TissueMultiscaleAnalysis::updateTime()
  {
    t = ((double)current_step+1.0)/num_load_steps;
  }
  int TissueMultiScaleAnalysis::run()
  {
    int result = 0;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    size_t cplng = getRelationID(amsi::getMultiscaleManager(),amsi::getScaleManager(),"macro","micro_fo");
#   ifdef LOGRUN
    initLogs();
#   endif
    std::vector<std::pair< double, double > > norm_hist;
    updateTime();
    tissue->setSimulationTime(t);
    logVolumes(vol_itms.begin(), vol_itms.end(), vols, current_step, tissue->getUField());
    tissue->computeInitGuess(las);
    tissue->initMicro();
    // Calculate constraints after initial guess
    tissue->updateConstraints();
    int complete = false;
    while(!complete)
    {
      std::cout << "Current load step = " << current_step << "." << std::endl;
      updateTime();
      tissue->setSimulationTime(t);
      if(amsi::numericalSolve(&itr,&cvg))
      {
        // needs to broadcast whether or not each iteration has converged
        cs->scaleBroadcast(cplng,&converged);
        if (current_step >= num_load_steps)
        {
          complete = true;
          std::cout << "Final load step converged! Simulation complete. Exiting..." << std::endl;
        }
        else
          tissue->step();
        current_step++;
      }
#     ifdef LOGRUN
      // displacement log
      logDisps(dsp_itms.begin(),dsp_itms.end(),disps,current_step,tissue->getUField());
      // force log
      for(auto mdl_ent = frc_itms.begin(); mdl_ent != frc_itms.end(); ++mdl_ent)
      {
        double frc[3] = {};
        tissue->getLoadOn((pGEntity)*mdl_ent,&frc[0]);
        if(rnk == 0)
          amsi::log(loads) << GEN_tag((pGEntity)*mdl_ent) << ", " << frc[0] << ", " << frc[1] << ", " << frc[2] << std::endl;
      }
      logVolumes(vol_itms.begin(), vol_itms.end(), vols, current_step, tissue->getUField());
      // write all streams
      std::fstream st_fs(state_file.c_str(), std::ios::out | std::ios::app);
      amsi::flushToStream(state,st_fs);
      if(rnk == 0)
      {
        std::fstream lds_fs(loads_file.c_str(), std::ios::out | std::ios::app);
        std::fstream dsp_fs(disps_file.c_str(), std::ios::out | std::ios::app);
        std::fstream vls_fs(vols_file.c_str(), std::ios::out | std::ios::app);
        std::fstream nrms_fs(norms_file.c_str(), std::ios::out | std::ios::app);
        std::fstream cnst_fs(cnstrnts_file.c_str(), std::ios::out | std::ios::app);
        amsi::flushToStream(loads,lds_fs);
        amsi::flushToStream(disps,dsp_fs);
        amsi::flushToStream(vols,vls_fs);
        amsi::flushToStream(norms,nrms_fs);
        amsi::flushToStream(cnstrnts,cnst_fs);
      }
#     endif
      // write mesh to file
      std::stringstream stpstrm;
      std::string pvd("/out.pvd");
      std::fstream pvdf;
      stpstrm << current_step;
      std::string msh_prfx("msh_stp_");
      apf::writeVtkFiles(std::string(amsi::fs->getResultsDir() + "/" + msh_prfx + stpstrm.str()).c_str(),tissue->getMesh());
      amsi::writePvdFile(pvd,msh_prfx,current_step);
      tissue->recoverSecondaryVariables(current_step);
      cs->scaleBroadcast(cplng,&complete);
    } // while(!complete)
    deleteLogs();
    return result;
  }
} // end of namespace Biotissue
