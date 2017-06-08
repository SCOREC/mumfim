#include "bioMultiscaleTissue.h"
#include "bioULMultiscaleIntegrator.h"
#include "bioULMultiscaleHydrostaticPressureIntegrator.h"
#include <bioRVEUtil.h> // micro
#include <amsiControlService.h> // amsi
namespace bio
{
  MultiscaleTissue::MultiscaleTissue(pGModel g, pParMesh m, pACase pd, MPI_Comm cm)
    : NonlinearTissue(g,m,pd,cm)
    , mltscl(NULL)
    , fo_cplg(apf_mesh,1) // get order from field
    , fbr_ornt(NULL)
    , nm_rves(0)
    , nm_rve_rgns(0)
    , rve_ents()
    , rve_ptrns()
    , snd_ptrns()
    , ini_ptrns()
    , rcv_ptrns()
  {
    fbr_ornt = apf::createIPField(apf_mesh,"P2",apf::SCALAR,1);
    mltscl = new ULMultiscaleIntegrator(&fo_cplg,apf_primary_field,crt_rve,1);
    M2m_id = amsi::getRelationID(amsi::getMultiscaleManager(),amsi::getScaleManager(),"macro","micro_fo");
    m2M_id = amsi::getRelationID(amsi::getMultiscaleManager(),amsi::getScaleManager(),"micro_fo","macro");
    apf::zeroField(fbr_ornt);
  }
  MultiscaleTissue::~MultiscaleTissue()
  {
    delete mltscl;
    apf::destroyField(fbr_ornt);
    apf::destroyField(prv_rve);
    apf::destroyField(crt_rve);
  }
  void MultiscaleTissue::Assemble(amsi::LAS * las)
  {
    computeRVEs();
    apfSimFEA::ApplyBC_Neumann(las);
    apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
    for(apf::MeshEntity * me = NULL; (me = apf_mesh->iterate(it));)
    {
      apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,me);
      amsi::ElementalSystem * sys = getIntegrator(me,0); // ERROR: assumes 1 type per ent
      sys->process(mlm);
      apf::NewArray<apf::Vector3> dofs;
      apf::getVectorNodes(sys->getElement(),dofs);
      apf::NewArray<int> ids;
      apf::getElementNumbers(apf_primary_numbering,me,ids);
      AssembleDOFs(las,
                   sys->numElementalDOFs(),
                   &ids[0],
                   &dofs[0],
                   &sys->getKe()(0,0),
                   &sys->getfe()(0),
                   sys->includesBodyForces());
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(it);
    for(auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->apply(las);
  }
  void MultiscaleTissue::computeRVEs()
  {
    std::vector<micro_fo_data> fo_data;
    serializeRVEData(std::back_inserter(fo_data)); // serialize data for fiber_only
    fo_cplg.sendRVEData(fo_data);
    fo_cplg.recvRVEData();
  }
  void MultiscaleTissue::initMicro()
  {
    /*
    amsi::ControlService * cs = amsi::ControlService::Instance();
    amsi::Task * macro = amsi::getLocal();
    amsi::DataDistribution * dd = amsi::createDataDistribution(macro,"micro_fo_data");
    (*dd) = 0;
    amsi::Assemble(dd,macro->comm());
    snd_ptrns[FIBER_ONLY] = cs->CreateCommPattern("micro_fo_data","macro","micro_fo");
    cs->CommPattern_Reconcile(snd_ptrns[FIBER_ONLY]);
    rcv_ptrns[FIBER_ONLY] = cs->RecvCommPattern("macro_fo_data","micro_fo","micro_fo_results","macro");
    cs->CommPattern_Reconcile(rcv_ptrns[FIBER_ONLY]);
    */
    computeRVETypeInfo();
    int num_rve_tps = rve_tps.size();
    cs->scaleBroadcast(M2m_id,&num_rve_tps);
    int rve_cnts[num_rve_tps];
    int ii = 0;
    for(auto tp = rve_tps.begin(); tp != rve_tps.end(); ++tp)
    {
      std::vector<MPI_Request> rqsts;
      cs->aSendBroadcast(std::back_inserter(rqsts),M2m_id,tp->first.c_str(),tp->first.size()+1);
      rve_cnts[ii++] = tp->second;
    }
    MPI_Request hdr_rqst;
    cs->aSendBroadcast(&hdr_rqst,M2m_id,&rve_cnts[0],num_rve_tps);
  }
  void MultiscaleTissue::updateMicro()
  {
    updateRVETypes();
    updateRVEExistence();
  }
  void MultiscaleTissue::updateRVETypes()
  {
    nm_rves = 0;
    apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
    for(apf::MeshEntity * me = NULL; (me = apf_mesh->iterate(it));)
    {
      apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,me);
      int ip = apf::countIntPoints(mlm,1);
      for(int ii = 0; ii < ip; ++ii)
      {
        apf::setScalar(prv_rve,me,ii,apf::getScalar(crt_rve,me,ii));
        int nw_tp = updateRVEType(me);
        apf::setScalar(crt_rve,me,ii,nw_tp);
        if(nw_tp == FIBER_ONLY)
          nm_rves++;
      }
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(it);
  }
  int MultiscaleTissue::updateRVEType(apf::MeshEntity * me)
  {
    // TODO: Use error estimate as determination of microscale TYPE
    //apf::MeshElement * me = apf::createMeshElement(apf_mesh,m_ent);
    //apf::Element * err_elmt = apf::createElement(apf_size_field,me);
    // get value of size field on element
    // if value is ?less? than some threshold, FIBER_ONLY, otherwise constitutive
    /*
      apf::NewArray<apf::Vector3> init_coords;
      apf::getVectorNodes(e_init_coords,init_coords);
      double x0 = 3.5;
      double pi = 4.0*atan(1.0);
      double x1 = 7.5 + 2*sin( 0.25*pi*(load_step+2) ) ;
      bool initial = (x0 < init_coords[0][0] && init_coords[0][0] < x1);
    */
    // Compute strain magnitude
    /*
      apf::Matrix3x3 eps;
      apf::getMatrix(strn,m_ent,0,eps);
      double strainMag = sqrt( pow(eps[0][0],2.0) +
      pow(eps[1][1],2.0) +
      pow(eps[2][2],2.0) +
      pow(eps[0][1],2.0) +
      pow(eps[1][2],2.0) +
      pow(eps[0][2],2.0) );
      if( false strainMag > 0.02 )
      {
      return FIBER_ONLY;
      }
      else
      {
      return NONE;
      }
    */
    apf::Matrix3x3 F;
    apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,me);
    apf::Element * e = apf::createElement(apf_primary_field,mlm);
    apf::Vector3 pcoords;
    int pt_num = 0; // assuming element with one integration point.
    apf::getIntPoint(mlm,1,pt_num,pcoords); // assume polynomial order of accuracy = 1.
    amsi::deformationGradient(e,pcoords,F);
    /*
    if (detF > 0.6)
      return FIBER_ONLY;
    else
      return NONE;
    */
    // multiscale based purely off of simmodeler specification, no adaptivity, need differentiate between initialization and updating
    pEntity snt = reinterpret_cast<pEntity>(me);
    pGEntity gsnt = EN_whatIn(snt);
    pAttribute mm = GEN_attrib(gsnt,"material model");
    pAttribute sm = Attribute_childByType(mm,"multiscale model");
    return  (sm ? FIBER_ONLY : NONE);
  }
  void MultiscaleTissue::updateRVEExistence()
  {
    std::vector<apf::MeshEntity*> nw_ents;
    std::vector<micro_fo_header> nw_hdrs;
    std::vector<micro_fo_params> nw_prms;
    std::vector<micro_fo_init_data> nw_data;
    std::vector<int> to_dlt;
    updateRVEDeletion(std::back_inserter(to_dlt));
    updateRVEAddition(std::back_inserter(nw_ents),
                      std::back_inserter(nw_hdrs),
                      std::back_inserter(nw_prms),
                      std::back_inserter(nw_data));
    std::copy(nw_ents.begin(),nw_ents.end(),std::back_inserter(rve_ents));
    // update comm patterns
    amsi::ControlService * cs = amsi::ControlService::Instance();
    amsi::Task * lt = amsi::getLocal();
    cs->RemoveData(snd_ptrns[FIBER_ONLY],to_dlt);
    std::vector<int> to_add(nw_ents.size());
    std::fill(to_add.begin(),to_add.end(),0);
    size_t add_id = cs->AddData(snd_ptrns[FIBER_ONLY],rve_ents,to_add);
    nm_rves += to_add.size();
    amsi::DataDistribution * dd = lt->getDD("micro_fo_data");
    (*dd) = nm_rves;
    amsi::Assemble(dd,lt->comm());
    cs->CommPattern_Assemble(snd_ptrns[FIBER_ONLY]);
    cs->Communicate(add_id,nw_hdrs,mtd.hdr);
    cs->Communicate(add_id,nw_prms,mtd.prm);
    cs->Communicate(add_id,nw_data,mtd.ini);
    /*
    // check for microscale load balancing (this reorders fiberOnly_elements depending on the microscale load balancing)
    cs->shareMigration(send_patterns[FIBER_ONLY],fiber_only_elements);
    */
    // update micro->macro communication
    cs->CommPattern_Reconcile(rcv_ptrns[FIBER_ONLY]);
  }
  amsi::ElementalSystem * MultiscaleTissue::getIntegrator(apf::MeshEntity * me, int ip)
  {
    int tp = apf::getScalar(crt_rve,me,ip);
    switch(tp)
    {
    case NONE:
      return constitutives[R_whatIn((pRegion)me)];
    case FIBER_ONLY:
      return mltscl;
    default:
      return NULL;
    }
  }
  void MultiscaleTissue::computeRVETypeInfo()
  {
    pGEntity rgn = NULL;
    GRIter ri = GM_regionIter(model);
    while((rgn = (pGEntity)GRIter_next(ri)))
    {
      pAttribute mdl = GEN_attrib(rgn,"material model");
      pAttribute sm = Attribute_childByType(mdl, "multiscale model");
      if(sm)
      {
        pAttributeString dir = (pAttributeString)Attribute_childByType(sm,"directory");
        pAttributeString prfx = (pAttributeString)Attribute_childByType(sm,"prefix");
        pAttributeInt cnt = (pAttributeInt)Attribute_childByType(sm,"count");
        char * dir_str = AttributeString_value(dir);
        char * tp_str = AttributeString_value(prfx);
        std::string tp(std::string(dir_str) + std::string("/") + std::string(tp_str));
        if(rve_tps.find(tp) == rve_tps.end())
          rve_tps[tp] = AttributeInt_value(cnt);
        Sim_deleteString(dir_str);
        Sim_deleteString(tp_str);
      }
    }
  }
  int MultiscaleTissue::getRVEType(apf::ModelEntity * ent)
  {
    int ii = -1;
    pGEntity rgn = reinterpret_cast<pGEntity>(ent);
    pAttribute mdl = GEN_attrib(rgn,"material model");
    pAttribute sm = Attribute_childByType(mdl, "multiscale model");
    if(sm)
    {
      pAttributeString dir = (pAttributeString)Attribute_childByType(sm,"directory");
      pAttributeString prfx = (pAttributeString)Attribute_childByType(sm,"prefix");
      char * dir_str = AttributeString_value(dir);
      char * tp_str = AttributeString_value(prfx);
      std::string tp(std::string(dir_str) + std::string("/") + std::string(tp_str));
      auto fnd = rve_tps.find(tp);
      if(fnd != rve_tps.end())
	ii = std::distance(rve_tps.begin(),fnd);
      Sim_deleteString(dir_str);
      Sim_deleteString(tp_str);
    }
    return ii;
  }
}
