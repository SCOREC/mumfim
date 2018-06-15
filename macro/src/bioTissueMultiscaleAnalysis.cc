#include "bioTissueMultiscaleAnalysis.h"
#include "bioAnalysisIO.h"
#include "bioMultiscaleTissue.h"
#include "bioMultiscaleConvergence.h"
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
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <fstream>
namespace bio
{
  void MultiscaleTissueIteration::iterate()
  {
    if(!PCU_Comm_Self())
      std::cout << "Multiscale Nonlinear Iteration : " << iteration() << std::endl;
    if(iteration() == 0)
      tssu->updateMicro();
    las->iter();
    fem_iter->iterate();
    tssu->iter();
    amsi::ModularIteration::iterate();
  }
  MultiscaleTissueAnalysis::MultiscaleTissueAnalysis(pGModel imdl, pParMesh imsh, pACase cs, MPI_Comm cm)
    : TissueAnalysis(imdl,imsh,cs,cm)
    , cplng(getRelationID(amsi::getMultiscaleManager(),amsi::getScaleManager(),"macro","micro_fo"))
  { }
  void MultiscaleTissueAnalysis::init()
  {
    // util data
    int rnk = -1;
    MPI_Comm_rank(cm,&rnk);
    // simmetrix attributes
    pACase pd = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
    pACase ss = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::SOLUTION_STRATEGY));
    // analysis params
    MultiscaleTissue * mt = NULL;
    tssu = mt = new MultiscaleTissue(mdl,msh,pd,cm);
    mx_stp = AttInfoInt_value((pAttInfoInt)AttNode_childByType((pANode)ss,"num timesteps"));
    dt = (double)1.0/(double)mx_stp;
    amsi::ModularIteration * mdl_itr = NULL;
    itr = mdl_itr = new MultiscaleTissueIteration(static_cast<MultiscaleTissue*>(tssu),las);
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
    cvg = new MultiscaleConvergence(cvg_stps.begin(),cvg_stps.end(),cplng);
    static_cast<MultiscaleTissue*>(tssu)->initMicro();
    // output params
#ifdef LOGRUN
    std::ostringstream cnvrt;
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
    amsi::log(state) << "STEP, ITER,   T, DESC" << std::endl;
#endif
  }
}
