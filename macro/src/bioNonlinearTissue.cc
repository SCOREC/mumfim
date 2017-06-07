#include "bioNonlinearTissue.h"
#include "bioAnalysis.h"
#include "bioNeoHookeanIntegrator.h"
#include "bioTrnsIsoNeoHookeanIntegrator.h"
#include "bioHolmesMowIntegrator.h"
//#include "StressStrainIntegrator.h"
#include "bioVariableRecovery.h"
#include <ErrorEstimators.h>
#include <apfFunctions.h>
#include <apfsimWrapper.h>
#include <amsiControlService.h>
#include <simBoundaryConditions.h>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <list>
#include <numeric>
#include <vector>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
namespace bio
{
  NonlinearTissue::NonlinearTissue(pGModel imdl, pParMesh imsh, pACase pd, MPI_Comm cm)
    : FEA(cm)
    , apfSimFEA(imdl,imsh,pd,cm)
    , constitutives()
    , dv_prev(0.0)
    , load_step(0)
    , iteration(0)
  {
    apf_primary_field = apf::createLagrangeField(apf_mesh,"displacement",apf::VECTOR,1);
    delta_u = apf::createLagrangeField(apf_mesh,"displacement_delta",apf::VECTOR,1);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);
    strs = apf::createIPField(apf_mesh,"stress",apf::MATRIX,1);
    rcvrd_strs = apf::createLagrangeField(apf_mesh,"recovered_stress",apf::MATRIX,1);
    strn = apf::createIPField(apf_mesh,"strain",apf::MATRIX,1);
    amsi::applyUniqueRegionTags(imdl,part,apf_mesh);
    std::vector<pANode> vol_cnst_nds;
    amsi::cutPaste<pANode>(AttNode_childrenByType((pANode)pd,"incompressible"),std::back_inserter(vol_cnst_nds));
    for(auto vol_nd = vol_cnst_nds.begin(); vol_nd != vol_cnst_nds.end(); ++vol_nd)
      vol_cnst.push_back(buildVolumeConstraint(pd,*vol_nd,apf_primary_numbering));
    pGEntity rgn = NULL;
    GRIter ri = GM_regionIter(imdl);
    while((rgn = (pGEntity)GRIter_next(ri)))
    {
      pAttribute mm = GEN_attrib(rgn,"material model");
      pAttribute cm = Attribute_childByType(mm,"continuum model");
      char * img_cls = Attribute_imageClass(cm);
      int cnst_type = getConstitutiveTypeFromString(img_cls);
      Sim_deleteString(img_cls);
      if(cnst_type == isotropic_neohookian) // push evaluation inside of elemental system constructor?
      {
        // should check to make sure the continuum model is iso lin ela for init solve?
        pAttributeTensor0 yng = (pAttributeTensor0)Attribute_childByType(cm,"youngs modulus");
        pAttributeTensor0 psn = (pAttributeTensor0)Attribute_childByType(cm,"poisson ratio");
        double E = AttributeTensor0_value(yng);
        double v = AttributeTensor0_value(psn);
        constitutives[rgn] = new NeoHookeanIntegrator(this,apf_primary_field,E,v,1);
      }
      else if(cnst_type == transverse_isotropic)
      {
        pAttributeTensor0 yng = (pAttributeTensor0)Attribute_childByType(cm,"youngs modulus");
        pAttributeTensor0 psn = (pAttributeTensor0)Attribute_childByType(cm,"poisson ratio");
        pAttributeTensor1 trnsvrs_axs = (pAttributeTensor1)Attribute_childByType(cm,"axis");
        pAttributeTensor0 trnsvrs_shr = (pAttributeTensor0)Attribute_childByType(cm,"axial shear modulus");
        pAttributeTensor0 trnsvrs_ygn = (pAttributeTensor0)Attribute_childByType(cm,"axial youngs modulus");
        double E = AttributeTensor0_value(yng);
        double v = AttributeTensor0_value(psn);
        double axs[3] = {0.0,0.0,0.0}; // axis must be constant for now
        AttributeTensor1_evalTensorDT(trnsvrs_axs,0.0,&axs[0]);
        double tG = AttributeTensor0_value(trnsvrs_shr);
        double tE = AttributeTensor0_value(trnsvrs_ygn);
        constitutives[rgn] = new TrnsIsoNeoHookeanIntegrator(this,apf_primary_field,E,v,&axs[0],tG,tE,1);
      }
    }
    GRIter_delete(ri);
    int dir_tps[] = {amsi::FieldUnit::displacement};
    amsi::buildSimBCQueries(pd,amsi::BCType::dirichlet,&dir_tps[0],(&dir_tps[0])+1,std::back_inserter(dir_bcs));
    int neu_tps[] = {amsi::NeuBCType::traction, amsi::NeuBCType::pressure};
    amsi::buildSimBCQueries(pd,amsi::BCType::neumann,&neu_tps[0],(&neu_tps[0])+2,std::back_inserter(neu_bcs));
  }
  NonlinearTissue::~NonlinearTissue()
  {
    apf::destroyField(delta_u);
    apf::destroyField(strs);
    apf::destroyField(rcvrd_strs);
    apf::destroyField(strn);
  }
  void NonlinearTissue::computeInitGuess(amsi::LAS * las)
  {
    LinearTissue lt(model,mesh,prob_def,analysis_comm);
    lt.setSimulationTime(T);
    LinearSolver(&lt,las);
    las->iter();
    apf::copyData(delta_u,lt.getField());
    apf::copyData(apf_primary_field,lt.getField());
  }
  void NonlinearTissue::ApplyBC_Dirichlet()
  {
    //apf::Field * nw_dlta_u = apf::createLagrangeField(apf_mesh,"tmp_delta_u",apf::VECTOR,1);
    //apf::copyData(nw_dlta_u,apf_primary_field);
    //amsi::PrintField(delta_u,std::cout).run();
    //amsi::PrintField(apf_primary_field,std::cout).run();
    // apply the new dirichlet bcs to the primary field
    amsi::apfSimFEA::ApplyBC_Dirichlet();
    //amsi::PrintField(apf_primary_field,std::cout).run();
    // subtract the current accumulated field (with new dirichlet bcs) from the old field to get new deltas
    //apf::axpy(-1.0,apf_primary_field,nw_dlta_u); // new - old = new_delta
    //amsi::PrintField(nw_dlta_u,std::cout).run();
    //amsi::WriteNZOp wrtr(1e-6);
    //amsi::MergeFields(delta_u,nw_dlta_u,&wrtr).run();
    //amsi::PrintField(delta_u,std::cout).run();
    //apf::destroyField(nw_dlta_u);
  }
  void NonlinearTissue::step()
  {
    for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->step();
    iteration = 0;
    load_step++;
  }
  void NonlinearTissue::iter()
  {
    for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->iter();
    iteration++;
  }
  void NonlinearTissue::Assemble(amsi::LAS * las)
  {
#   ifdef LOGRUN
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
    amsi::Log macro_efficiency = amsi::activateLog("macro_efficiency");
    double pre_assmbl = amsi::getElapsedTime(macro_efficiency);
    amsi::log(macro_efficiency)  << load_step << ", " << iteration << ", " << pre_assmbl << ", ACTIVE, PRE_ASSEMBLE" << std::endl;
#   endif
    ApplyBC_Neumann(las);
    apf::MeshEntity * me = NULL;
    // custom iterator would be perfect for switching for multiscale version
    auto it = apf_mesh->begin(analysis_dim);
    while((me = apf_mesh->iterate(it)))
    {
      amsi::ElementalSystem * constitutive = constitutives[R_whatIn((pRegion)me)];
      apf::MeshElement * mlmt = apf::createMeshElement(apf_mesh,me);
      //amsi::ElementalSystem * sys = getElementalSystem(me,0); // assumes 1 type of system per element
      constitutive->process(mlmt);
      apf::NewArray<apf::Vector3> dofs;
      apf::getVectorNodes(constitutive->getElement(),dofs);
      apf::NewArray<int> ids;
      apf::getElementNumbers(apf_primary_numbering,me,ids);
      AssembleDOFs(las,
                   constitutive->numElementalDOFs(),
                   &ids[0],
                   &dofs[0],
                   &constitutive->getKe()(0,0),
                   &constitutive->getfe()(0),
                   constitutive->includesBodyForces());
      apf::destroyMeshElement(mlmt);
    }
    apf_mesh->end(it);
    double nrm = 0.0;
    las->GetVectorNorm(nrm);
    // process constraints
    for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->apply(las);
    las->GetVectorNorm(nrm);
#   ifdef LOGRUN
    double post_assmbl = amsi::getElapsedTime(macro_efficiency);
    amsi::log(macro_efficiency)  << load_step << ", " << iteration << ", " << post_assmbl << ", IDLE, POST_ASSEMBLE" << std::endl;
    double pre_slv = amsi::getElapsedTime(macro_efficiency);
    amsi::log(macro_efficiency)  << load_step << ", " << iteration << ", " << pre_slv << ", ACTIVE, PRE_SOLVE" << std::endl;
#   endif
  }
  void NonlinearTissue::UpdateDOFs(const double * sol)
  {
#   ifdef LOGRUN
    amsi::Log macro_efficiency = amsi::activateLog("macro_efficiency");
    double post_slv = amsi::getElapsedTime(macro_efficiency);
    amsi::log(macro_efficiency)  << load_step << ", " << iteration << ", " << post_slv << ", ACTIVE, POST_SOLVE" << std::endl;
#   endif
    // accumulate displacement deltas into primary field
    //apfSimFEA::UpdateDOFs(sol);
    amsi::AccumOp ac_op;
    amsi::FreeApplyOp frac_op(apf_primary_numbering,&ac_op);
    amsi::ApplyVector(apf_primary_numbering,apf_primary_field,sol,first_local_dof,&frac_op).run();
    //amsi::PrintField(apf_primary_field,std::cout).run();
    amsi::WriteOp wr_op;
    amsi::FreeApplyOp frwr_op(apf_primary_numbering,&wr_op);
    amsi::ApplyVector(apf_primary_numbering,delta_u,sol,first_local_dof,&frwr_op).run();
    //amsi::PrintField(delta_u,std::cout).run();
    apf::synchronize(apf_primary_field);
    apf::synchronize(delta_u);
  }
  // this is just a field op... should be replacable with what we currently provide in amsi
  void NonlinearTissue::computeDispL2Norm(double & norm)
  {
    norm = 0.0;
    double sqrtnorm;
    // iterate over all mesh entities, from vertices up.. (if needed)
    int field_components = apf::countComponents(apf_primary_field);
    for(int ii = 0; ii < 3; ii++)
    {
      apf::FieldShape * fs = apf::getShape(apf_primary_field);
      if(!fs->hasNodesIn(ii))
        break;
      apf::MeshIterator * it = apf_mesh->begin(ii);
      while(apf::MeshEntity * me = apf_mesh->iterate(it))
      {
        apf::Vector3 incr_disp;
        int num_nodes = fs->countNodesOn(apf_mesh->getType(me));
        for(int jj = 0 ; jj < num_nodes; jj++)
        {
          apf::getVector(delta_u,me,jj,incr_disp);
          for(int kk = 0; kk < field_components; kk++)
          {
            // Norm is sum of incremental displacements
            sqrtnorm = incr_disp[kk];
            norm += sqrtnorm * sqrtnorm;
          }
        }
      }
    }
    // Compute norm over all local processes
    norm = sqrt(amsi::comm_sum(norm));
  }
  void NonlinearTissue::getLoadOn(pGEntity ent, double * frc)
  {
    amsi::SimBCQuery * snbcq = amsi::findSimBCQueryOn(neu_bcs.begin(),neu_bcs.end(),ent);
    assert(snbcq);
    assert(!snbcq->isSpaceExpr());
    for(int ii = 0; ii < 3; ii++)
      frc[ii] = snbcq->getValue(ii,T);
  }
  void NonlinearTissue::recoverSecondaryVariables(int load_step)
  {
    //#ifdef SCOREC
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
    std::stringstream fnm;
    fnm << amsi::fs->getResultsDir()
        << "/qlty.stp_" << load_step << ".rnk_" << rnk << ".dat";
    // analyze and print the quality of the elements
    apf::Field * qfld = amsi::analyzeMeshQuality(apf_mesh,apf_primary_field);
    std::ofstream file(fnm.str().c_str(),std::ofstream::out);
    amsi::PrintField(qfld,file).run();
    apf::destroyField(qfld);
    //#endif
  }
  void NonlinearTissue::storeStrain(apf::MeshElement * me, double* strain)
  {
    apf::MeshEntity * m_ent = apf::getMeshEntity(me);
    apf::Matrix3x3 eps(strain[0],strain[3],strain[5],
                       strain[3],strain[1],strain[4],
                       strain[5],strain[4],strain[2]);
    apf::setMatrix(strn,m_ent,0,eps);
  }
  void NonlinearTissue::storeStress(apf::MeshElement * me, double* stress)
  {
    apf::MeshEntity * m_ent = apf::getMeshEntity(me);
    apf::Matrix3x3 sigma(stress[0],stress[3],stress[5],
                         stress[3],stress[1],stress[4],
                         stress[5],stress[4],stress[2]);
    apf::setMatrix(strs,m_ent,0,sigma);
  }
}
