#include "bioMultiscaleTissue.h"
#include "bioULMultiscaleIntegrator.h"
#include "bioULMultiscaleHydrostaticPressureIntegrator.h"
#include <amsiControlService.h> // amsi
#include <amsiDetectOscillation.h>
namespace bio
{
  MultiscaleTissue::MultiscaleTissue(pGModel g, pParMesh m, pACase pd,
                                     pACase ss, MPI_Comm cm)
      : NonlinearTissue(g, m, pd, ss, cm)
      , mltscl(NULL)
      , crt_rve(apf::createIPField(apf_mesh, "micro_rve_type", apf::SCALAR, 1))
      , prv_rve(apf::createIPField(apf_mesh, "micro_old_type", apf::SCALAR, 1))
      , compute_ornt_3D(false)
      , compute_ornt_2D(false)
      , ornt_3D(NULL)
      , ornt_2D(NULL)
      , fo_cplg(apf_mesh, crt_rve, prv_rve,
                apf::getShape(apf_primary_field)->getOrder())
      , nm_rves(0)
      , rve_dirs()
  {
    ornt_3D = apf::createIPField(apf_mesh, "ornt_tens_3D", apf::MATRIX,1);
    ornt_2D = apf::createIPField(apf_mesh, "ornt_tens_2D", apf::MATRIX,1);
    mltscl = new ULMultiscaleIntegrator(&fo_cplg,strn,strs,apf_primary_field,dfm_grd,1);
    M2m_id = amsi::getRelationID(amsi::getMultiscaleManager(),amsi::getScaleManager(),"macro","micro_fo");
    m2M_id = amsi::getRelationID(amsi::getMultiscaleManager(),amsi::getScaleManager(),"micro_fo","macro");
    apf::zeroField(ornt_3D);
    apf::zeroField(ornt_2D);
    apf::zeroField(crt_rve);
    apf::zeroField(prv_rve);
    // get properties for output orienation tensor
    pAttribute ornt_tens = GM_attrib(g, "output orientation tensor");
    pAttribute ornt_tens_3D = NULL;
    pAttribute ornt_tens_2D = NULL;
    if(ornt_tens)
    {
      ornt_tens_3D =
          Attribute_childByType(ornt_tens, "3D Orientation Tensor");
      ornt_tens_2D =
          Attribute_childByType(ornt_tens, "2D Orientation Tensor");
    }
    // set the class orientation status
    compute_ornt_3D = ornt_tens_3D != NULL;
    compute_ornt_2D = ornt_tens_2D != NULL;
    pAttributeTensor1 axs = NULL;
    if (compute_ornt_2D) {
      axs = (pAttributeTensor1)Attribute_childByType(ornt_tens_2D, "axis");
    }
    ornt_2D_axis[0] = compute_ornt_2D ? AttributeTensor1_value(axs, 0) : 0.0;
    ornt_2D_axis[1] = compute_ornt_2D ? AttributeTensor1_value(axs, 1) : 0.0;
    ornt_2D_axis[2] = compute_ornt_2D ? AttributeTensor1_value(axs, 2) : 0.0;
  }
  MultiscaleTissue::~MultiscaleTissue()
  {
    delete mltscl;
    apf::destroyField(crt_rve);
    apf::destroyField(prv_rve);
    apf::destroyField(ornt_3D);
    apf::destroyField(ornt_2D);
  }
  void MultiscaleTissue::Assemble(amsi::LAS * las)
  {
    computeRVEs();
#ifdef LOGRUN
    amsi::Log state = amsi::activateLog("tissue_efficiency");
    amsi::log(state) << load_step << ", "
                     << iteration << ", "
                     << MPI_Wtime() << ", "
                     << "start_fea"
                     << std::endl;
#endif
    apfSimFEA::ApplyBC_Neumann(las);
    apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
    for(apf::MeshEntity * me = NULL; (me = apf_mesh->iterate(it));)
    {
      apf::MeshElement * mlm = apf::createMeshElement(current_coords,me);
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
#ifdef LOGRUN
    amsi::log(state) << load_step << ", "
                     << iteration << ", "
                     << MPI_Wtime() << ", "
                     << "end_fea"
                     << std::endl;
    amsi::log(state) << load_step << ", "
                     << iteration << ", "
                     << MPI_Wtime() << ", "
                     << "start_solve"
                     << std::endl;
#endif
  }
  void MultiscaleTissue::computeRVEs()
  {
#ifdef LOGRUN
    amsi::Log state = amsi::activateLog("tissue_efficiency");
    amsi::log(state) << load_step << ", "
                     << iteration << ", "
                     << MPI_Wtime() << ", "
                     << "start_rves"
                     << std::endl;
#endif

    std::vector<micro_fo_data> fo_data;
    serializeRVEData(std::back_inserter(fo_data)); // serialize data for fiber_only
    fo_cplg.sendRVEData(fo_data);
    fo_cplg.recvRVEData();
#ifdef LOGRUN
    amsi::log(state) << load_step << ", "
                     << iteration << ", "
                     << MPI_Wtime() << ", "
                     << "end_rves"
                     << std::endl;
#endif
  }
  void MultiscaleTissue::initMicro()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    fo_cplg.initCoupling();
    loadRVELibraryInfo();
    int num_rve_tps = rve_dirs.size();
    cs->scaleBroadcast(M2m_id,&num_rve_tps);
    std::vector<int> rve_cnts(num_rve_tps);
    int ii = 0;
    for(auto tp = rve_dirs.begin(); tp != rve_dirs.end(); ++tp)
    {
      std::vector<MPI_Request> rqsts;
      cs->aSendBroadcast(std::back_inserter(rqsts),M2m_id,tp->c_str(),tp->size()+1);
      rve_cnts[ii++] = rve_dir_cnts[std::distance(rve_dirs.begin(),tp)];
    }
    std::vector<MPI_Request> hdr_rqst;
    cs->aSendBroadcast(std::back_inserter(hdr_rqst),M2m_id,&rve_cnts[0],num_rve_tps);
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
    apf::destroyElement(e);
    apf::destroyMeshElement(mlm);
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
    std::vector<micro_fo_header> nw_hdrs;
    std::vector<micro_fo_params> nw_prms;
    std::vector<micro_fo_init_data> nw_data;
    std::vector<int> to_dlt;
    std::vector<micro_fo_solver> slvr_prms;
    std::vector<micro_fo_int_solver> slvr_int_prms;
    fo_cplg.updateRVEDeletion(std::back_inserter(to_dlt));
    serializeNewRVEData(std::back_inserter(nw_hdrs),
                        std::back_inserter(nw_prms),
                        std::back_inserter(nw_data),
                        std::back_inserter(slvr_prms),
                        std::back_inserter(slvr_int_prms));
    fo_cplg.deleteRVEs(to_dlt.begin(), to_dlt.end());
    size_t add_ptrn = fo_cplg.addRVEs(nw_hdrs.size());
    fo_cplg.sendNewRVEs(add_ptrn, nw_hdrs, nw_prms, nw_data, slvr_prms,
                        slvr_int_prms);
    fo_cplg.updateRecv();
    nm_rves += nw_hdrs.size();
  }
  void MultiscaleTissue::recoverSecondaryVariables(int step)
  {
    NonlinearTissue::recoverSecondaryVariables(step);
    fo_cplg.recvRVEStepData();
    // loop over all mesh regions
    apf::MeshIterator* it = apf_mesh->begin(analysis_dim);
    while (apf::MeshEntity* meshEnt = apf_mesh->iterate(it)) {
      apf::MeshElement* mlm = apf::createMeshElement(apf_mesh, meshEnt);
      int numIP = apf::countIntPoints(mlm, 1);
      for (int ip = 0; ip < numIP; ++ip) {
        micro_fo_step_result* stp_rslt = fo_cplg.getRVEStepResult(meshEnt, 0);
        // this is kindof ugly due to repetitions, but it takes out some
        // logic from the tight loop
        if (compute_ornt_2D && compute_ornt_3D) {
          apf::Matrix3x3 orn_tens_3D;
          apf::Matrix3x3 orn_tens_2D;
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              orn_tens_3D[i][j] = stp_rslt->data[i * 3 + j];
              orn_tens_2D[i][j] = stp_rslt->data[9 + i * 3 + j];
            }
          }
          apf::setMatrix(ornt_3D, meshEnt, ip, orn_tens_3D);
          apf::setMatrix(ornt_2D, meshEnt, ip, orn_tens_2D);
        }
        else if (compute_ornt_3D) {
          apf::Matrix3x3 orn_tens_3D;
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              orn_tens_3D[i][j] = stp_rslt->data[i * 3 + j];
            }
          }
          apf::setMatrix(ornt_3D, meshEnt, ip, orn_tens_3D);
        }
        else if (compute_ornt_2D) {
          apf::Matrix3x3 orn_tens_2D;
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              orn_tens_2D[i][j] = stp_rslt->data[9 + i * 3 + j];
            }
          }
          apf::setMatrix(ornt_2D, meshEnt, ip, orn_tens_2D);
        }
      }
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(it);
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
  void MultiscaleTissue::loadRVELibraryInfo()
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
        if(std::find(rve_dirs.begin(), rve_dirs.end(), tp) == rve_dirs.end())
        {
          rve_dirs.push_back(tp);
          rve_dir_cnts.push_back(AttributeInt_value(cnt));
        }
        Sim_deleteString(dir_str);
        Sim_deleteString(tp_str);
      }
    }
  }
  void MultiscaleTissue::getExternalRVEData(apf::MeshEntity* ent,
                                            micro_fo_header& hdr,
                                            micro_fo_params& prm,
                                            micro_fo_solver& slvr,
                                            micro_fo_int_solver& int_slvr)
  {
    pGEntity smdl_ent = EN_whatIn(reinterpret_cast<pEntity>(ent));
    pAttribute mm = GEN_attrib(smdl_ent, "material model");
    pAttribute sm = Attribute_childByType(mm, "multiscale model");
    assert(sm);
    pAttributeTensor0 fbr_rd =
        (pAttributeTensor0)Attribute_childByType(sm, "radius");
    pAttributeTensor0 vl_frc =
        (pAttributeTensor0)Attribute_childByType(sm, "volume fraction");
    pAttribute prms = Attribute_childByType(sm, "force reaction");
    pAttributeTensor0 yngs =
        (pAttributeTensor0)Attribute_childByType(prms, "youngs modulus");
    pAttributeTensor0 nnlr = (pAttributeTensor0)Attribute_childByType(
        prms, "nonlinearity parameter");
    pAttributeTensor0 lntr =
        (pAttributeTensor0)Attribute_childByType(prms, "linear transition");
    // get properties for output orienation tensor
    hdr.data[COMPUTE_ORIENTATION_3D] = compute_ornt_3D;
    hdr.data[COMPUTE_ORIENTATION_2D] = compute_ornt_2D;
    prm.data[ORIENTATION_AXIS_X] = ornt_2D_axis[0];
    prm.data[ORIENTATION_AXIS_Y] = ornt_2D_axis[1];
    prm.data[ORIENTATION_AXIS_Z] = ornt_2D_axis[2];
    prm.data[FIBER_RADIUS] = AttributeTensor0_value(fbr_rd);
    prm.data[VOLUME_FRACTION] = AttributeTensor0_value(vl_frc);
    prm.data[YOUNGS_MODULUS] = AttributeTensor0_value(yngs);
    prm.data[NONLINEAR_PARAM] = nnlr ? AttributeTensor0_value(nnlr) : 0.0;
    prm.data[LINEAR_TRANSITION] = lntr ? AttributeTensor0_value(lntr) : 0.0;
    pAttribute micro_cnvg = AttCase_attrib(solution_strategy, "microscale convergence operator");
    assert(micro_cnvg);
    pAttribute dtct_osc =
        Attribute_childByType(micro_cnvg, "oscillation detection");
    assert(dtct_osc);
    pAttribute micro_cnvg_tol = Attribute_childByType(micro_cnvg, "micro convergence tolerance");
    assert(micro_cnvg_tol);
    pAttribute micro_slvr_tol = Attribute_childByType(micro_cnvg, "micro solver tolerance");
    pAttribute num_attempts = NULL;
    pAttribute cut_factor = NULL;
    pAttribute itr_cap = NULL;
    pAttribute prev_itr_factor = NULL;
    char* detect_osc_type = Attribute_imageClass(dtct_osc);
    if (strcmp(detect_osc_type, "iteration only") == 0) {
      num_attempts = Attribute_childByType(dtct_osc, "number of attempts");
      cut_factor = Attribute_childByType(dtct_osc, "attempt cut factor");
      itr_cap = Attribute_childByType(dtct_osc, "micro iteration cap");
      int_slvr.data[DETECT_OSCILLATION_TYPE] =
          static_cast<int>(amsi::DetectOscillationType::IterationOnly);
      assert(itr_cap);
    }
    else if (strcmp(detect_osc_type, "previous residual") == 0) {
      num_attempts = Attribute_childByType(dtct_osc, "number of attempts");
      cut_factor = Attribute_childByType(dtct_osc, "attempt cut factor");
      prev_itr_factor = Attribute_childByType(dtct_osc, "previous itr factor");
      int_slvr.data[DETECT_OSCILLATION_TYPE] =
          static_cast<int>(amsi::DetectOscillationType::PrevNorm);
      assert(prev_itr_factor);
    }
    else if (strcmp(detect_osc_type, "combined") == 0) {
      num_attempts = Attribute_childByType(dtct_osc, "number of attempts");
      cut_factor = Attribute_childByType(dtct_osc, "attempt cut factor");
      itr_cap = Attribute_childByType(dtct_osc, "micro iteration cap");
      prev_itr_factor = Attribute_childByType(dtct_osc, "previous itr factor");
      int_slvr.data[DETECT_OSCILLATION_TYPE] =
          static_cast<int>(amsi::DetectOscillationType::IterationPrevNorm);
      assert(itr_cap);
      assert(prev_itr_factor);
    }
    else {
      std::cerr << detect_osc_type
                << " is not a valid oscillation detection type. attDefs and "
                   "MultiscaleTissue are out of sync.\n";
      std::abort();
    }
    delete [] detect_osc_type;
    detect_osc_type = NULL;
    assert(num_attempts);
    assert(cut_factor);
    // FIXME replace before commit
    slvr.data[MICRO_CONVERGENCE_TOL] =
        AttributeDouble_value((pAttributeDouble)micro_cnvg_tol);
    slvr.data[MICRO_SOLVER_TOL] =
        micro_slvr_tol ? AttributeDouble_value((pAttributeDouble)micro_slvr_tol)
                       : 1E-6;
    // should not choose number less than 1E-6 for now due to LAS using single
    // precision
    assert(slvr.data[MICRO_SOLVER_TOL] <= slvr.data[MICRO_CONVERGENCE_TOL]);
    slvr.data[PREV_ITER_FACTOR] =
        prev_itr_factor
            ? AttributeDouble_value((pAttributeDouble)prev_itr_factor)
            : 0;
    int_slvr.data[MAX_MICRO_CUT_ATTEMPT] =
        AttributeInt_value((pAttributeInt)num_attempts);
    int_slvr.data[MICRO_ATTEMPT_CUT_FACTOR] =
        AttributeInt_value((pAttributeInt)cut_factor);
    int_slvr.data[MAX_MICRO_ITERS] =
        itr_cap ? AttributeInt_value((pAttributeInt)itr_cap) : 0;
  }
  void MultiscaleTissue::getInternalRVEData(apf::MeshEntity * rgn,
                                            micro_fo_header & hdr,
                                            micro_fo_params & ,
                                            micro_fo_init_data & dat)
  {
    apf::ModelEntity * mnt = reinterpret_cast<apf::ModelEntity*>(EN_whatIn(reinterpret_cast<pEntity>(rgn)));
    hdr.data[RVE_TYPE]     = getRVEDirectoryIndex(rve_dirs.begin(),rve_dirs.end(),mnt);
    hdr.data[FIELD_ORDER]  = apf::getShape(delta_u)->getOrder();
    hdr.data[ELEMENT_TYPE] = apf_mesh->getType(rgn);
    hdr.data[GAUSS_ID]     = -1; // needs to be set for each IP
    apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,rgn);
    apf::Element * ce = apf::createElement(apf_mesh->getCoordinateField(),mlm);
    apf::NewArray<apf::Vector3> Ni;
    apf::getVectorNodes(ce,Ni);
    int nds = apf::countNodes(ce);
    for(int nd = 0; nd < nds; ++nd)
      Ni[nd].toArray(&dat.init_data[nd*3]);
    apf::destroyElement(ce);
    apf::destroyMeshElement(mlm);
  }
}
