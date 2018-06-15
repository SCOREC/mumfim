#include "bioTissueAnalysis.h"
#include "bioAnalysisIO.h"
namespace bio
{
  TissueAnalysis::TissueAnalysis(pGModel md, pParMesh ms, pACase ca, MPI_Comm c)
    : cm(c)
    , mdl(md)
    , msh(ms)
    , cs(ca)
    , t(0.0)
    , dt(0.0)
    , stp(0)
    , mx_stp(1)
    , tssu(NULL)
    , itr()
    , cvg()
    , cvg_stps()
    , trkd_vols()
    , las(new amsi::PetscLAS(0,0))
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
  { }
  TissueAnalysis::~TissueAnalysis()
  {
    delete cvg;
    for(auto c = cvg_stps.begin(); c != cvg_stps.end(); ++c)
    {
      amsi::Convergence * cn = *c;
      delete cn;
      *c = NULL;
    }
    delete tssu;
    delete las;
  }
  void TissueAnalysis::init()
  {
    // util data
    int rnk = -1;
    MPI_Comm_rank(cm,&rnk);
    // simmetrix attributes
    pACase pd = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
    pACase ss = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::SOLUTION_STRATEGY));
    // analysis params
    tssu = new NonlinearTissue(mdl,msh,pd,cm);
    amsi::ModularIteration * mdl_itr = NULL;
    itr = mdl_itr = new TissueIteration(tssu,las);
    mx_stp = AttInfoInt_value((pAttInfoInt)AttNode_childByType((pANode)ss,"num timesteps"));
    dt = (double)1.0/(double)mx_stp;
    std::vector<pANode> trk_vols;
    amsi::cutPaste<pANode>(AttNode_childrenByType((pANode)ss,"track volume"),std::back_inserter(trk_vols));
    std::vector<VolCalc*> vls;
    for(auto trk_vol = trk_vols.begin(); trk_vol != trk_vols.end(); ++trk_vol)
    {
      std::vector<apf::ModelEntity*> mdl_ents;
      amsi::getAssociatedModelItems(ss,*trk_vol,std::back_inserter(mdl_ents));
      trkd_vols[*trk_vol] = new VolCalc(mdl_ents.begin(),mdl_ents.end(),tssu->getUField());
      mdl_itr->addOperation(trkd_vols[*trk_vol]);
    }
    buildLASConvergenceOperators(ss,itr,las,std::back_inserter(cvg_stps));
    buildVolConvergenceOperators(ss,itr,tssu->getUField(),trkd_vols,std::back_inserter(cvg_stps));
    cvg_stps.push_back(new amsi::ResetIteration(&amsi::linear_convergence, itr));
    cvg = new amsi::MultiConvergence(cvg_stps.begin(),cvg_stps.end());
    // output params
#ifdef LOGRUN
    std::stringstream cnvrt;
    cnvrt << rnk;
    state_fn = amsi::fs->getResultsDir() + "/tissue_state." + cnvrt.str() + ".log";
    amsi::getTrackedModelItems(cs,"output force",std::back_inserter(frc_itms));
    amsi::getTrackedModelItems(cs,"output displacement",std::back_inserter(dsp_itms));
    amsi::getTrackedModelItems(cs,"output volume",std::back_inserter(vol_itms));
    // initialize logging
    state = amsi::activateLog("tissue_efficiency");
    if(rnk == 0)
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
    t += dt;
    tssu->setSimulationTime(t);
    logVolumes(vol_itms.begin(),vol_itms.end(), vols, stp, tssu->getUField());
    tssu->computeInitGuess(las);
    completed = false;
    while(!completed)
    {
#ifdef LOGRUN
      amsi::log(state) << stp << ", "
                       << itr->iteration() << ", "
                       << MPI_Wtime() << ", "
                       << "start_step"
                       << std::endl;
#endif
      if(!PCU_Comm_Self())
	std::cout << "Load step = " << stp << std::endl;
      if(amsi::numericalSolve(itr,cvg))
      {
#ifdef LOGRUN
        amsi::log(state) << stp << ", "
                         << itr->iteration() << ", "
                         << MPI_Wtime() << ", "
                         << "end_solve"
                         << std::endl;
#endif
        if(stp == mx_stp)
        {
          completed = true;
          std::cout << "Final load step converged. Case complete." << std::endl;
        }
        else
        {
          for(auto vol = trkd_vols.begin(); vol != trkd_vols.end(); ++vol)
            vol->second->step();
          las->step();
          tssu->step();
        }
#ifdef LOGRUN
      amsi::log(state) << stp << ", "
                       << itr->iteration() << ", "
                       << MPI_Wtime() << ", "
                       << "end_step"
                       << std::endl;
#endif
      }
      else
      {
        completed = true;
        std::cerr << "ERROR: Step " << stp << " failed to converge!" << std::endl;
      }
      logDisps(dsp_itms.begin(),dsp_itms.end(),dsps,stp,tssu->getUField());
      logForces(frc_itms.begin(),frc_itms.end(),frcs,stp,tssu);
      logVolumes(vol_itms.begin(),vol_itms.end(),vols,stp,tssu->getUField());
      stp++;
      t += dt;
      tssu->setSimulationTime(t);
      checkpoint();
    }
    deinit();
  }
  void TissueAnalysis::checkpoint()
  {
    int rnk = -1;
    MPI_Comm_rank(cm,&rnk);
    if(rnk == 0)
    {
      std::ofstream frcs_fs(frcs_fn.c_str(), std::ios::out | std::ios::app);
      std::ofstream dsps_fs(dsps_fn.c_str(), std::ios::out | std::ios::app);
      std::ofstream vols_fs(vols_fn.c_str(), std::ios::out | std::ios::app);
      std::ofstream nrms_fs(nrms_fn.c_str(), std::ios::out | std::ios::app);
      std::ofstream cnst_fs(constraint_fn.c_str(), std::ios::out | std::ios::app);
      amsi::flushToStream(frcs,frcs_fs);
      amsi::flushToStream(dsps,dsps_fs);
      amsi::flushToStream(vols,vols_fs);
      amsi::flushToStream(nrms,nrms_fs);
      amsi::flushToStream(constraints,cnst_fs);
    }
    std::ofstream st_fs(state_fn.c_str(), std::ios::out | std::ios::app);
    amsi::flushToStream(state,st_fs);
    // write mesh to file
    std::string pvd("/out.pvd");
    std::ofstream pvdf;
    std::string msh_prfx("msh_stp_");
    std::stringstream cnvrt;
    cnvrt << stp;
    apf::writeVtkFiles(std::string(amsi::fs->getResultsDir() + "/" + msh_prfx + cnvrt.str()).c_str(),tssu->getMesh());
    amsi::writePvdFile(pvd,msh_prfx,stp);
    tssu->recoverSecondaryVariables(stp);
  }
  void TissueAnalysis::revert()
  {

  }
  void TissueAnalysis::deinit()
  {
    int rnk = -1;
    MPI_Comm_rank(cm,&rnk);
    if(rnk == 0)
    {
      amsi::deleteLog(vols);
      amsi::deleteLog(dsps);
      amsi::deleteLog(nrms);
      amsi::deleteLog(frcs);
      amsi::deleteLog(constraints);
    }
    amsi::deleteLog(state);
  }
}
