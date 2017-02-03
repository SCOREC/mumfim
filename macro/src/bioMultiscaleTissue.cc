#include "bioMultiscaleTissue.h"
#include "bioULMultiscaleIntegrator.h"
#include <RVE_Util.h> // micro
#include <amsiControlService.h> // amsi
namespace bio
{
  MultiscaleTissue::MultiscaleTissue(pGModel g, pParMesh m, pACase pd, MPI_Comm cm)
    : NonlinearTissue(g,m,pd,cm)
    , mltscl(NULL)
    , crt_rve(NULL)
    , prv_rve(NULL)
    , fbr_ornt(NULL)
    , nm_rves(0)
    , nm_rve_rgns(0)
    , rve_ents()
    , rslts()
    , rslt_mp()
    , rve_ptrns()
    , snd_ptrns()
    , ini_ptrns()
    , rcv_ptrns()
    , mtd()
  {
    // primary field created in NonlinearTissue
    crt_rve = apf::createIPField(apf_mesh,"current_rve",apf::SCALAR,1);
    prv_rve = apf::createIPField(apf_mesh,"previous_rve",apf::SCALAR,1);
    fbr_ornt = apf::createIPField(apf_mesh,"fiber_orientation",apf::MATRIX,1);
    mltscl = new ULMultiscaleIntegrator(this,apf_primary_field,crt_rve,1);
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
      (*cnst)->apply(las,apf_mesh,part,apf_primary_numbering);
  }
  void MultiscaleTissue::computeRVEs()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    std::vector<micro_fo_data> fo_data;
    rslts.clear(); // make output param
    rslt_mp.clear();
    serializeRVEData(std::back_inserter(fo_data)); // serialize data for fiber_only
    if(lb_per_iteration) // on first iteration (0) this should simply be empty...
    {
      cs->shareMigration(snd_ptrns[FIBER_ONLY],rve_ents);
      cs->CommPattern_Reconcile(rcv_ptrns[FIBER_ONLY]);
    }
    cs->Communicate(snd_ptrns[FIBER_ONLY],fo_data,mtd.micro_fo_data_type);    // send fiber_only data
    cs->Communicate(rcv_ptrns[FIBER_ONLY],rslts,mtd.micro_fo_result_type);    // wait for fiber_only results
    std::list<apf::MeshEntity*>::iterator me = rve_ents.begin(); // apply result data (implicit ordering)
    int ent_rve = 0;
    int ii = 0;
    for(std::vector<micro_fo_result>::iterator it = rslts.begin(); it != rslts.end(); ++it)
    {
      int ent_rves = countRVEsOn(*me);
      std::vector<RVE_Info> & rslt_vc = rslt_mp[*me]; // retrieve the result vector related to the current element
      RVE_Info & info = rslt_vc[ent_rve];             // set result struct info
      info.mType = apf_mesh->getType(*me);
      info.order = 1;
      info.derivS = &it->data[0];
      // todo (m) : fix hacky hard-coded bs
      apf::setMatrix(fbr_ornt,*me,0,apf::Matrix3x3(it->data[81],it->data[82],it->data[83],
                                                   it->data[84],it->data[85],it->data[86],
                                                   it->data[87],it->data[88],it->data[89]));
      if(ent_rve == ent_rves-1)
      {
        ent_rve = 0;
        me++;
      }
      else
      {
        ent_rve++;
      }
      ii++;
    }
  }
  void MultiscaleTissue::initMicro()
  {
    mtd.MultiscaleDataTypesMPICommit();
    amsi::ControlService * cs = amsi::ControlService::Instance();
    amsi::Task * macro = amsi::getLocal();
    macro->createDD("micro_fo_data");
    macro->setLocalDDValue("micro_fo_data",0);
    macro->assembleDD("micro_fo_data");
    snd_ptrns[FIBER_ONLY] = cs->CreateCommPattern("micro_fo_data","macro","micro_fo");
    cs->CommPattern_Reconcile(snd_ptrns[FIBER_ONLY]);
    rcv_ptrns[FIBER_ONLY] = cs->RecvCommPattern("macro_fo_data","micro_fo","micro_fo_results","macro");
    cs->CommPattern_Reconcile(rcv_ptrns[FIBER_ONLY]);
  }
  int MultiscaleTissue::countRVEsOn(apf::MeshEntity * me)
  {
    int cnt = 0;
    if(apf::getScalar(crt_rve,me,0) == FIBER_ONLY)
    {
      apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,me);
      cnt = apf::countIntPoints(mlm,getOrder(mlm));
      apf::destroyMeshElement(mlm);
    }
    return cnt;
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
    lt->setLocalDDValue("micro_fo_data",nm_rves);
    lt->assembleDD("micro_fo_data");
    cs->CommPattern_Assemble(snd_ptrns[FIBER_ONLY]);
    cs->Communicate(add_id,nw_hdrs,mtd.micro_fo_header_data_type);
    cs->Communicate(add_id,nw_prms,mtd.micro_fo_parameter_data_type);
    cs->Communicate(add_id,nw_data,mtd.micro_fo_init_data_type);
    /*
    // check for microscale load balancing (this reorders fiberOnly_elements depending on the microscale load balancing)
    cs->shareMigration(send_patterns[FIBER_ONLY],fiber_only_elements);
    */
    // update micro->macro communication
    cs->CommPattern_Reconcile(rcv_ptrns[FIBER_ONLY]);
  }
  RVE_Info * MultiscaleTissue::getRVEResult(apf::MeshEntity * me, int ip)
  {
    return &rslt_mp[me][ip]; // check existence
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
}
