#include "bioNonlinearTissue.h"
#include "bioBatheMultiscale.h"
#include "bioMultiscaleIntegrator.h"
#include "bioULMultiscaleIntegrator.h"
#include "bioNeoHookeanIntegrator.h"
#include "bioHolmesMowIntegrator.h"
#include "StressStrainIntegrator.h"
#include "bioPreserveVolConstraintIntegrator.h"
#include "RVE_Util.h"
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
    , poisson_ratio(0.236249)
    , youngs_modulus(43200000)
    , dv_prev(0.0)
    , load_step(0)
    , iteration(0)
    , iteration_beta(0)
  {
    apf_primary_field = apf::createLagrangeField(apf_mesh,"displacement",apf::VECTOR,1);
    delta_u = apf::createLagrangeField(apf_mesh,"displacement_delta",apf::VECTOR,1);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);
    strs = apf::createIPField(apf_mesh,"stress",apf::MATRIX,1);
    rcvrd_strs = apf::createLagrangeField(apf_mesh,"recovered_stress",apf::MATRIX,1);
    strn = apf::createIPField(apf_mesh,"strain",apf::MATRIX,1);
    amsi::applyUniqueRegionTags(imdl,part,apf_mesh);
    pGEntity rgn = NULL;
    GRIter ri = GM_regionIter(imdl);
    while((rgn = (pGEntity)GRIter_next(ri)))
    {
      pAttribute inc = GEN_attrib(rgn,"specify incompressible");
      if(inc)
      {
        // which version (surface, volume, etc)
        // enforcement method
        // parameters
        double beta = 0.0;
        pAttribute mthd = Attribute_childByType(inc,"incompressible enforcement");
        char * mtd = Attribute_imageClass(mthd);
        if(std::string("penalty method").compare(mtd) == 0)
        {
          pAttributeDouble bt_att = (pAttributeDouble)Attribute_childByType(mthd,"beta");
          beta = AttributeDouble_value(bt_att);
        }
        Sim_deleteString(mtd);
        pAttributeInt vrsn_att = (pAttributeInt)Attribute_childByType(inc,"version");
        if(AttributeInt_value(vrsn_att) == 0)
        {
          pPList fcs = GR_faces((pGRegion)rgn);
          void * itr = 0;
          pGFace fc;
          while((fc = (pGFace)PList_next(fcs,&itr)))
          {
            auto cnst = new VolumeConstraintSurface(fc,GEN_tag((pGRegion)rgn),part,delta_u,1);
            cnst->setBeta(beta);
            vol_cnst.push_back(cnst);
          }
          PList_delete(fcs);
        }
        else
          std::cout << "WARNING: unsupported incompressibility version specified!" << std::endl;
      }
      pAttribute mm = GEN_attrib(rgn,"material model");
      pAttribute cm = Attribute_childByType(mm,"continuum model");
      // should check to make sure the continuum model is iso lin ela for init solve?
      pAttributeTensor0 yngs = (pAttributeTensor0)Attribute_childByType(cm,"youngs modulus");
      pAttributeTensor0 psn = (pAttributeTensor0)Attribute_childByType(cm,"poisson ratio");
      youngs_modulus = AttributeTensor0_value(yngs);
      poisson_ratio = AttributeTensor0_value(psn);
      stress_strain_system = new amsi::LinearStressStrainIntegrator(apf_primary_field,
                                                                    strn,
                                                                    strs,
                                                                    youngs_modulus,poisson_ratio);
    }
    GRIter_delete(ri);
    double shear_modulus = ( 3.0 * youngs_modulus * (1.0 - 2.0 * poisson_ratio) )/( 2.0 * (1.0 + poisson_ratio) );
    constitutive =  new NeoHookeanIntegrator(this,apf_primary_field,shear_modulus,poisson_ratio,1);
    int dir_tps[] = {amsi::FieldUnit::displacement};
    amsi::buildSimBCQueries(pd,amsi::BCType::dirichlet,&dir_tps[0],(&dir_tps[0])+1,std::back_inserter(dir_bcs));
    int neu_tps[] = {amsi::NeuBCType::traction, amsi::NeuBCType::pressure};
    amsi::buildSimBCQueries(pd,amsi::BCType::neumann,&neu_tps[0],(&neu_tps[0])+2,std::back_inserter(neu_bcs));
    /// Initialize rgn_vols vectors.
    calcInitVolumes(); ///< init_rgn_vols initialized here.
    double vol_temp = 0.0;
    for (uint ii=0; ii<init_rgn_vols.size(); ii++)
      vol_temp += init_rgn_vols[ii];
    // Initialize initial and current (global) volumes for each instant of VolumeConstraintSurface class.
    for(std::vector<VolumeConstraintSurface*>::iterator cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++){
      (*cnst)->setInitVol(vol_temp);
      (*cnst)->setVol(vol_temp);
      (*cnst)->setPrevVol(vol_temp);
    }
    rgn_vols = init_rgn_vols;
    prev_rgn_vols = init_rgn_vols;
  }
  NonlinearTissue::~NonlinearTissue()
  {
    apf::destroyField(delta_u);
    apf::destroyField(strs);
    apf::destroyField(rcvrd_strs);
    apf::destroyField(strn);
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
  void NonlinearTissue::RenumberDOFs()
  {
    //constraint_dofs = vol_cnst.size(); // implicit input the renumbering function
    amsi::apfFEA::RenumberDOFs();
    //int cntr = 1;
    /*
      for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->setDof(global_dof_count - cntr); // implicit output from the renumbering function
    */
  }
  void NonlinearTissue::step()
  {
    //apf::zeroField(delta_u);
    //amsi::displaceMesh(delta_u);
    //amsi::displaceMesh(apf_primary_field);
    iteration = 0;
    iteration_beta = 0;
    // maybe update the neumann bcs values HERE?
    load_step++;
    //apf::zeroField(apf_primary_field);
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
    // process constraints
    for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->apply(las,apf_mesh,part,apf_primary_numbering);

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
  void NonlinearTissue::ComputeDispL2Norm(double & norm)
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
    std::stringstream filename;
    filename << amsi::fs->getResultsDir()
             << "mesh_quality_step" << load_step << "." << rnk << ".dat";
    // analyze and print the quality of the elements
    apf::Field * quality_field = amsi::analyzeMeshQuality(apf_mesh,apf_primary_field);
    std::ofstream file(filename.str().c_str(),std::ofstream::out);
    amsi::PrintField(quality_field,file).run();
    apf::destroyField(quality_field);
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
  void NonlinearTissue::updateConstraints()
  {
    // Update beta and lambda based on values of dv and v calculated above.
    double max_beta = 128;
    double max_iter = 15;
    double idx = 0;
    double eps = 1e-2;
    for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
    {
      double dv = rgn_vols[idx] - init_rgn_vols[idx]; // dv is calculated with respect to initial (target) volume.
      double dv_prev = prev_rgn_vols[idx] - init_rgn_vols[idx];
      if (std::abs(dv) > eps * init_rgn_vols[idx])    // criteria for volume difference constraint.
      {
        // beta update
        if (std::abs(dv_prev) > 1e-10 && dv > 0.5 * dv_prev)
        {
          double beta = 0.0;
          if (iteration < max_iter && (*cnst)->getBeta() < max_beta)
          {
            beta = (*cnst)->getBeta() * 2.0;
            (*cnst)->setBeta(beta);
          }
          else
          {
            beta = (*cnst)->getBeta() * 0.5;
            (*cnst)->setBeta(beta);
          }
        }
        // Update lambda
        std::cout << "Beta = " << (*cnst)->getBeta() << std::endl;
        //(*cnst)->setLambda((*cnst)->getLambda() + (dv * (*cnst)->getBeta()));
        (*cnst)->setLambda(0.0);
        std::cout << "Lambda = " << (*cnst)->getLambda() << std::endl;
      }
      idx++;
    }
  }
  void NonlinearTissue::updateConstraintsAccm()
  {
    // Determine accumulated volume change
    double dv = 0.0;
    double v = 0.0;
    double v0 = 0.0;
    int rgns = numVolumeConstraints();
    for (int ii=0; ii<rgns; ii++)
    {
      v += rgn_vols[ii];
      v0 += init_rgn_vols[ii];
    }
    dv = v - v0;
    // Update beta and lambda based on values of dv and v calculated above.
    for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
    {
      // Update lambda
      std::cout << "Beta = " << (*cnst)->getBeta() << std::endl;
      //(*cnst)->setLambda((*cnst)->getLambda() + (dv/v0 * (*cnst)->getBeta()));
      (*cnst)->setLambda(0.0);
      std::cout << "dv = " << dv << std::endl;
      std::cout << "Lambda = " << (*cnst)->getLambda() << std::endl;
      (*cnst)->setGflag(true);
    }
  }
  void NonlinearTissue::updateConstraintsAccm_Incrmt()
  {
    // Determine accumulated volume change
    double dv = 0.0;
    double vp = 0.0;
    double dv_rgn = 0;
    int rgns = numVolumeConstraints();
    for (int ii=0; ii<rgns; ii++)
    {
      dv_rgn = rgn_vols[ii] - prev_rgn_vols[ii];
      dv += dv_rgn;
      vp += prev_rgn_vols[ii];
    }
    // Update beta and lambda based on values of dv and v calculated above.
    for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
    {
      // Update lambda
      //std::cout << "Beta = " << (*cnst)->getBeta() << std::endl;
      (*cnst)->setLambda((*cnst)->getLambda() + (dv/vp * (*cnst)->getBeta()));
      //(*cnst)->setLambda(0.0);
      //std::cout << "dv = " << dv << std::endl;
      //std::cout << "Lambda = " << (*cnst)->getLambda() << std::endl;
    }
    dv_prev = dv;
  }
  void NonlinearTissue::logCnstrntParams(int ldstp, int iteration, int rnk)
  {
    amsi::Log cnstrnts = amsi::activateLog("constraints");
    double lmbda = 0.0;
    double beta = 0.0;
    if (rnk == 0 && vol_cnst.size() > 0){
      lmbda = vol_cnst[0]->getLambda();
      beta = vol_cnst[0]->getBeta();
      amsi::log(cnstrnts) << ldstp << ", "
                          << iteration << ", "
                          << lmbda << ", "
                          << beta << std::endl;
    }
  }
} // end of namespace Biotissue
