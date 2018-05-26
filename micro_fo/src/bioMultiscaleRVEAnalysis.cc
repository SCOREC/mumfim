#include "bioMultiscaleRVEAnalysis.h"
#include "bioFiberRVEAnalysis.h"
#include "bioFiberNetworkIO.h"
#include "bioMicroFOMultiscale.h"
#include "bioRVEVolumeTerms.h"
#include <apfFEA.h> // amsi
#include <apfMatrixUtil.h> //amsi
#include <apfMeshIterator.h> // amsi
#include <apfMDS.h>
#include <gmi.h>
#include <cassert>
namespace bio
{
  FiberRVEAnalysis * initFromMultiscale(FiberNetwork * fn,
                                        las::CSR * csr,
                                        las::SparskitBuffers * bfrs,
                                        micro_fo_header & hdr,
                                        micro_fo_params & prm,
                                        micro_fo_init_data & ini)
  {
    FiberRVEAnalysis * rve = makeFiberRVEAnalysis(fn,csr,bfrs);
    rve->multi = new MultiscaleRVE(rve->rve,fn,hdr,prm,ini);
    return rve;
  }
  class ApplyDeformationGradient : public amsi::FieldOp
  {
  protected:
    apf::Field * xyz;
    apf::Field * du;
    apf::Field * u;
    apf::MeshEntity * ent;
    apf::Matrix3x3 FmI;
  public:
    ApplyDeformationGradient(apf::Matrix3x3 F, apf::Mesh * msh, apf::Field * du_, apf::Field * u_)
      : xyz(msh->getCoordinateField())
      , du(du_)
      , u(u_)
      , ent(NULL)
    {
      int d = msh->getDimension();
      for(int ii = 0; ii < d; ++ii)
        for(int jj = 0; jj < d; ++jj)
          FmI[ii][jj] = F[ii][jj] - (ii == jj ? 1.0 : 0.0);
    }
    virtual bool inEntity(apf::MeshEntity * m)
    {
      ent = m;
      return true;
    }
    virtual void outEntity() {}
    virtual void atNode(int nd)
    {
      apf::Vector3 nd_xyz;
      apf::Vector3 nd_u_old;
      apf::getVector(xyz,ent,nd,nd_xyz);
      apf::getVector(u,ent,nd,nd_u_old);
      apf::Vector3 nd_u = FmI * nd_xyz;
      apf::Vector3 nd_du = nd_u - nd_u_old;
      apf::setVector(u,ent,nd,nd_u);
      apf::setVector(du,ent,nd,nd_du);
    }
    void run()
    {
      apply(u);
    }
  };
  void applyMultiscaleCoupling(FiberRVEAnalysis * ans, micro_fo_data * data)
  {
    int d = ans->rve->getDim();
    assert(d == 3 || d == 2);
    apf::Matrix3x3 F;
    for(int ei = 0; ei < d; ++ei)
      for(int ej = 0; ej < d; ++ej)
        F[ei][ej] = data->data[ei*d + ej];
    apf::Mesh * rve_msh = ans->rve->getMesh();
    apf::Field * rve_du = ans->rve->getdUField();
    apf::Field * rve_u = ans->rve->getUField();
    ApplyDeformationGradient(F,rve_msh,rve_du,rve_u).run();
    apf::Mesh * fn_msh = ans->fn->getNetworkMesh();
    apf::Field * fn_du = ans->fn->getdUField();
    apf::Field * fn_u = ans->fn->getUField();
    ApplyDeformationGradient(F,fn_msh,fn_du,fn_u).run();
  }
  // might be duplicated from bioFiberRVEAnalysis.cc ?
  void recoverMicroscaleStress(FiberRVEAnalysis * ans, double * stress)
  {
    int dim = ans->fn->getNetworkMesh()->getDimension();
    apf::Field * xyz = ans->fn->getNetworkMesh()->getCoordinateField();
    apf::Vector3 xyz_val;
    apf::Vector3 f_val;
    apf::Numbering * num = ans->fn->getUNumbering();
    int dofs[3] {};
    double * f = NULL;
    ans->ops->get(ans->f,f);
    apf::Matrix3x3 strs(0.0,0.0,0.0,
                        0.0,0.0,0.0,
                        0.0,0.0,0.0);
    for(auto nd = ans->bnd_nds.begin(); nd != ans->bnd_nds.end(); ++nd)
    {
      apf::getVector(xyz,*nd,0,xyz_val);
      for(int ii = 0; ii < ans->fn->getDim(); ++ii)
        dofs[ii] = apf::getNumber(num,*nd,0,ii);
      strs[0][0] += xyz_val[0] * f[dofs[0]];
      strs[0][1] += xyz_val[1] * f[dofs[0]];
      strs[0][2] += xyz_val[2] * f[dofs[0]];
      strs[1][0] += xyz_val[0] * f[dofs[1]];
      strs[1][1] += xyz_val[1] * f[dofs[1]];
      strs[1][2] += xyz_val[2] * f[dofs[1]];
      strs[2][0] += xyz_val[0] * f[dofs[2]];
      strs[2][1] += xyz_val[1] * f[dofs[2]];
      strs[2][2] += xyz_val[2] * f[dofs[2]];
    }
    // take the symmetric part of the current stress matrix
    apf::Matrix3x3 sym_strs = amsi::symmetricPart(strs);
    amsi::mat2VoigtVec(dim,sym_strs,&stress[0]);
  }
  void convertStress(FiberRVEAnalysis * ans, double * stress)
  {
    int dim = ans->fn->getNetworkMesh()->getDimension();
    int sigma_length = dim == 3 ? 6 : 3;
    // convert to a macro-scale term
    double vol = ans->rve->measureDu();
    double scale_conversion = ans->multi->getScaleConversion();
    for(int ii = 0; ii < sigma_length; ++ii)
      stress[ii] *= scale_conversion / vol;
  }
  void recoverAvgVolStress(FiberRVEAnalysis * ans, double * Q)
  {
    // Q terms
    (void)ans;
    (void)Q;
  }
  // honestly break most of this out into an class (maybe the drve/dfe one?) and find a way to require that the operations happen in the correct order (modifications of the stiffness matrix)
  void recoverStressDerivs(FiberRVEAnalysis * ans, double * sigma, double * dstrss_drve)
  {
    (void)ans;
    (void)dstrss_drve;
    // stress derivative terms
    apf::Mesh * fn_msh = ans->fn->getNetworkMesh();
    int dim = ans->rve->getMesh()->getDimension();
    int fn_dof_cnt = ans->fn->getDofCount();
    int rve_dof_cnt = ans->rve->numNodes()*dim;
    // wierd area term matrix
    // this is calculated on the reference domain and makes assumptions based on that fact
    apf::DynamicMatrix dRdx_rve(fn_dof_cnt,rve_dof_cnt);
    dRdx_rve.zero();
    apf::Numbering * rve_dofs = ans->rve->getNumbering();
    apf::Numbering * fn_dofs = ans->fn->getUNumbering();
    for(auto vrt = ans->bnd_nds.begin(); vrt != ans->bnd_nds.end(); ++vrt)
    {
      std::vector<apf::MeshEntity*> bnds[RVE::side::all];
      for(int sd = RVE::side::bot; sd != RVE::side::all; ++sd)
      {
        int plnr_dim = (sd == RVE::side::rgt || sd == RVE::side::lft) ? 0 : (sd == RVE::side::bot || sd == RVE::side::top) ? 1 : 2;
        RVE::side rve_sd = static_cast<RVE::side>(sd);
        apf::Vector3 crd;
        fn_msh->getPoint(*vrt,0,crd);
        if(ans->rve->onBoundary(crd,rve_sd))
        {
          apf::MeshEntity * sd_ent = ans->rve->getSide(rve_sd);
          double a = apf::measure(ans->rve->getMesh(),sd_ent);
          // extract method
          int vrt_cnt = dim == 3 ? 4 : 2;
          apf::MeshEntity * fc_vrts[vrt_cnt];
          ans->rve->getMesh()->getDownward(sd_ent,0,&fc_vrts[0]);
          apf::NewArray<int> rve_dof_ids(vrt_cnt * dim);
          apf::getElementNumbers(rve_dofs,sd_ent,rve_dof_ids);
          apf::NewArray<int> fn_dof_ids(dim);
          apf::getElementNumbers(fn_dofs,*vrt,fn_dof_ids);
          // iterate over the line/face (dep on dim)
          // assumes dim = 3
          // get the vert opposite the current vert on the face, which is what the
          //  value given by the calculation below is applicable to,
          //  we are basically calculating barycentric coordinates using
          //  the fact the RVE is unit cube as a huge shortcut inverting the global mapping
          int opp_vrt[] = {2,3,0,1};
          for(int ci = 0; ci < vrt_cnt; ++ci)
          {
            apf::Vector3 ci_crd;
            ans->rve->getMesh()->getPoint(fc_vrts[ci],0,ci_crd);
            // this method of calculating the area/length associated with the vertex is only valid on the reference domain as it depends on the verts of a face/edge being both coplanar and axis-aligned
            apf::Vector3 spn = ci_crd - crd;
            double ai = 1.0;
            for(int ii = 0; ii < dim; ++ii)
              ai *= ii == plnr_dim ? 1.0 : spn[ii];
            ai = fabs(ai / a);
            // iterate over the dim
            for(int ej = 0; ej < dim; ++ej)
              dRdx_rve(fn_dof_ids[ej],rve_dof_ids[opp_vrt[ci]*dim + ej]) = ai;
          }
          bnds[sd].push_back(*vrt);
          break;
        }
      }
    }
    // have dRdx_rve
    // setup additional solves with the rows of dRdx_rve as the force vector
    // need to modify matrix as done in calc_precond in old code prior to these solves, this also produces dS_dx_fn the change of the stresses on the boundary of the RVE w.r.t. the fiber network coordinates
    int sigma_length = dim == 3 ? 6 : 3;
    apf::DynamicMatrix dS_dx_fn(sigma_length,fn_dof_cnt);
    dS_dx_fn.zero();
    apf::Numbering * dofs = ans->fn->getUNumbering();
    double * F = NULL;
    ans->ops->get(ans->f,F);
    // iterate over all fibers with a node on the boundary
    for(auto vrt = ans->bnd_nds.begin(); vrt != ans->bnd_nds.end(); ++vrt)
    {
      apf::Mesh * fn_msh = ans->fn->getNetworkMesh();
      assert(fn_msh->countUpward(*vrt) == 1);
      apf::Vector3 Nx;
      // need the displaced node coordinates, not the reference coords
      apf::getVector(ans->fn->getXpUField(),*vrt,0,Nx);
      //fn_msh->getPoint(*vrt,0,Nx);
      apf::MeshEntity * edg = fn_msh->getUpward(*vrt,0);
      apf::MeshEntity * vrts[2];
      fn_msh->getDownward(edg,0,&vrts[0]);
      int bnd_dofs[dim];
      for(int dd = 0; dd < dim; ++dd)
        bnd_dofs[dd] = apf::getNumber(dofs,*vrt,0,dd);
      for(int vrt = 0; vrt < 2; ++vrt)
      {
        for(int dd = 0; dd < dim; ++dd)
        {
          int dof = apf::getNumber(dofs,vrts[vrt],0,dd);
          apf::Vector3 ks;
          ks.zero();
          for(int ii = 0; ii < dim; ++ii)
          {
            ks[ii] = las::getSparskitMatValue(ans->k,bnd_dofs[ii],dof);
            ks[ii] *= -1.0;
          }
          apf::Vector3 delta;
          for(int ii = 0; ii < dim; ++ii)
            delta[ii] = dof == bnd_dofs[ii] ? 1.0 : 0.0;
          apf::Vector3 fs;
          for(int ii = 0; ii < dim; ++ii)
            fs[ii] = F[bnd_dofs[ii]];
          apf::Matrix3x3 m = apf::tensorProduct(ks,Nx) + apf::tensorProduct(delta,fs);
          apf::Matrix3x3 ms = amsi::symmetricPart(m);
          apf::DynamicVector vls(sigma_length);
          amsi::mat2VoigtVec(dim,ms,&vls(0));
          for(int ii = 0; ii < sigma_length; ++ii)
            dS_dx_fn(ii,dof) += vls(ii);
          for(int ii = 0; ii < dim; ++ii)
            las::setSparskitMatValue(ans->k,bnd_dofs[ii],dof,0.0);
        }
      }
      for(int ii = 0; ii < dim; ++ii)
        las::setSparskitMatValue(ans->k,bnd_dofs[ii],bnd_dofs[ii],1.0);
    }
    // have dS_dx_fn
    apf::DynamicVector f(fn_dof_cnt);
    apf::DynamicVector u(fn_dof_cnt);
    las::Vec * skt_f = las::createSparskitVector(fn_dof_cnt);
    las::Vec * skt_u = las::createSparskitVector(fn_dof_cnt);
    las::LasSolve * slv = las::createSparskitQuickLUSolve(ans->slv);
    // zero rows in the matrix w/ boundary conditions
    // this effects the force vector which is used in the calculation of
    // dS_dx_fn above so we must do it after.
    // this might not be necessary anymore since the change in the
    //  matrix modification indices in the dS_dx_fn term
    applyRVEBC(ans->bnd_nds.begin(),ans->bnd_nds.end(),
               ans->fn->getUNumbering(),ans->ops,ans->k,ans->f);
    apf::DynamicMatrix dx_fn_dx_rve(fn_dof_cnt,rve_dof_cnt);
    for(int ii = 0; ii < rve_dof_cnt; ++ii)
    {
      // apf -> double * -> sparskit
      dRdx_rve.getColumn(ii,f);
      double * fptr = NULL;
      ans->ops->get(skt_f,fptr);
      std::copy(f.begin(),f.end(),fptr);
      // solve k u = f for modified f
      slv->solve(ans->k,skt_u,skt_f);
      // sparskit -> double * -> apf
      double * uptr = NULL;
      ans->ops->get(skt_u,uptr);
      std::copy(uptr,uptr+fn_dof_cnt,u.begin());
      dx_fn_dx_rve.setColumn(ii,u);
    }
    // have dx_fn_dx_rve
    apf::DynamicMatrix dS_dx_rve;
    apf::DynamicMatrix odS_dx_fn(sigma_length,fn_dof_cnt);
    apf::DynamicVector rw(fn_dof_cnt);
    int sigma_pmt[6] = {0, 3, 5, 4, 2, 1};
    for(int ii = 0; ii < sigma_length; ++ii)
    {
      dS_dx_fn.getRow(ii,rw);
      odS_dx_fn.setRow(sigma_pmt[ii],rw);
    }
    apf::DynamicMatrix odx_fn_dx_rve(fn_dof_cnt,rve_dof_cnt);
    apf::DynamicVector cl(fn_dof_cnt);
    int rve_pmt[8] = {2, 0, 6, 3, 4, 1, 7, 5};
    for(int ii = 0; ii < rve_dof_cnt; ++ii)
    {
      dx_fn_dx_rve.getColumn(ii,cl);
      odx_fn_dx_rve.setColumn((rve_pmt[ii/3]*3)+(ii%3),cl);
    }
    apf::multiply(odS_dx_fn,odx_fn_dx_rve,dS_dx_rve);
    apf::multiply(dS_dx_fn,dx_fn_dx_rve,dS_dx_rve);
    //apf::createMeshElement(ans->rve->getMesh(),ans->rve->getMeshEnt());
    apf::DynamicVector dV_dx_rve;
    double vol = ans->rve->measureDu();
    CalcdV_dx_rve calcdv_dx_rve(2,ans->rve->getUField());
    apf::MeshElement * mlm = apf::createMeshElement(ans->rve->getXpUField(),ans->rve->getMeshEnt());
    calcdv_dx_rve.process(mlm);
    calcdv_dx_rve.getdVdxrve(dV_dx_rve);
    apf::destroyMeshElement(mlm);
    double scale_conversion = ans->multi->getScaleConversion();
    apf::DynamicVector col(sigma_length);
    apf::Vector<6> vsig(sigma);
    for(int ii = 0; ii < rve_dof_cnt; ++ii)
    {
      apf::DynamicVector dsig = apf::fromVector(vsig);
      dS_dx_rve.getColumn(ii,col);
      col /= vol;
      dsig *= (dV_dx_rve[ii] / (vol * vol));
      col -= dsig;
      col *= scale_conversion;
      dS_dx_rve.setColumn(ii,col);
    }
    apf::DynamicMatrix dx_rve_dx_fe;
    ans->multi->calcdRVEdFE(dx_rve_dx_fe,ans->rve);
    apf::DynamicMatrix dS_dx;
    apf::multiply(dS_dx_rve,dx_rve_dx_fe,dS_dx);
    amsi::mat2Array(dS_dx,dstrss_drve);
  }
  void recoverMultiscaleResults(FiberRVEAnalysis * ans, micro_fo_result * data)
  {
    // rebuild everything since we want the force vector without
    // boundary conditions anyway
    ans->ops->zero(ans->k);
    apf::Mesh * fn  = ans->fn->getNetworkMesh();
    apf::MeshEntity * me = NULL;
    apf::MeshIterator * it = fn->begin(1);
    while((me = fn->iterate(it)))
    {
      apf::MeshElement * mlm = apf::createMeshElement(fn,me);
      ans->es->process(mlm);
      apf::destroyMeshElement(mlm);
    }
    fn->end(it);
    double * S = &data->data[0];
    recoverMicroscaleStress(ans,S);
    double * Q = &data->data[6];
    recoverAvgVolStress(ans,Q);
    double * dS_dx_rve = &data->data[9];
    // calculate macroscale stress derivs
    // requires microscale stress
    recoverStressDerivs(ans,S,dS_dx_rve);
    // convert microscale stress to macroscale
    convertStress(ans,S);
  }
  MultiscaleRVEAnalysis::MultiscaleRVEAnalysis()
    : eff()
    , wgt()
    , tmg()
    , recv_ptrn()
    , send_ptrn()
    , rve_dd(NULL)
    , M2m_id()
    , m2M_id()
//    , dat_tp()
    , rve_tp_cnt(0)
    , fns()
    , sprs()
    , bfrs(NULL)
    , macro_iter(0)
    , macro_step(0)
  {
    M2m_id = amsi::getRelationID(amsi::getMultiscaleManager(),amsi::getScaleManager(),"macro","micro_fo");
    m2M_id = amsi::getRelationID(amsi::getMultiscaleManager(),amsi::getScaleManager(),"micro_fo","macro");
  }
  void MultiscaleRVEAnalysis::init()
  {
    initLogging();
    initCoupling();
    initAnalysis();
  }
  void MultiscaleRVEAnalysis::initLogging()
  {
    eff = amsi::activateLog("micro_fo_efficiency");
    wgt = amsi::activateLog("micro_fo_weights");
    tmg = amsi::activateLog("micro_fo_timing");
  }
  void MultiscaleRVEAnalysis::initCoupling()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    rve_dd = amsi::createDataDistribution(amsi::getLocal(),"macro_fo_data");
    recv_ptrn = cs->RecvCommPattern("micro_fo_data","macro",
                                    "macro_fo_data","micro_fo");
    cs->CommPattern_Reconcile(recv_ptrn);
    send_ptrn = cs->CommPattern_UseInverted(recv_ptrn,"macro_fo_data",
                                            "micro_fo","macro");
    cs->CommPattern_Assemble(send_ptrn);
    cs->CommPattern_Reconcile(send_ptrn);
  }
  void MultiscaleRVEAnalysis::initAnalysis()
  {
    int num_rve_tps = 0;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->scaleBroadcast(M2m_id,&num_rve_tps);
    char ** rve_tp_dirs = new char * [num_rve_tps];
    MPI_Request rqsts[num_rve_tps];
    // the order of receipt might be non-deterministic. need to handle that
    for(int ii = 0; ii < num_rve_tps; ++ii)
    {
      int cnt = 0;
      while((cnt = cs->aRecvBroadcastSize<char>(M2m_id)) == 0)
      { }
      rve_tp_dirs[ii] = new char [cnt];
      // don't have to block to wait since we know the message was available for size info
      cs->aRecvBroadcast(&rqsts[ii],M2m_id,&rve_tp_dirs[ii][0],cnt);
    }
    MPI_Status stss[num_rve_tps];
    MPI_Waitall(num_rve_tps,&rqsts[0],&stss[0]);
    MPI_Request hdr_rqst;
    rve_tp_cnt.resize(num_rve_tps);
    cs->aRecvBroadcast(&hdr_rqst,M2m_id,&rve_tp_cnt[0],num_rve_tps);
    MPI_Status hdr_sts;
    MPI_Waitall(1,&hdr_rqst,&hdr_sts);
    // Read in all the fiber network meshes and reactions
    for(int ii = 0; ii < num_rve_tps; ++ii)
    {
      fns.push_back(new FiberNetworkReactions* [rve_tp_cnt[ii]]);
      sprs.push_back(new las::CSR* [rve_tp_cnt[ii]]);
    }
    int dof_max = -1;
    for(int ii = 0; ii < num_rve_tps; ii++)
    {
      for(int jj = 0; jj < rve_tp_cnt[ii]; ++jj)
      {
        std::stringstream fl;
        fl << rve_tp_dirs[ii] << jj+1 << ".txt";
        FiberNetworkReactions * fn_rctns = new FiberNetworkReactions;
        fn_rctns->msh = loadFromFile(fl.str());
        apf::Mesh2 * fn = fn_rctns->msh;
        fl << ".params";
        loadParamsFromFile(fn,fl.str(),std::back_inserter(fn_rctns->rctns));
        fns[ii][jj] = fn_rctns;
        apf::Field * u = apf::createLagrangeField(fn,"u",apf::VECTOR,1);
        apf::zeroField(u);
        apf::Numbering * n = apf::createNumbering(u);
        int dofs = apf::NaiveOrder(n);
        sprs[ii][jj] = las::createCSR(n,dofs);
        apf::destroyNumbering(n);
        apf::destroyField(u);
        dof_max = dofs > dof_max ? dofs : dof_max;
      }
    }
    bfrs = new las::SparskitBuffers(dof_max);
  }
  void MultiscaleRVEAnalysis::updateCoupling()
  {
    std::vector<micro_fo_header> hdrs;
    std::vector<micro_fo_params> prms;
    std::vector<micro_fo_init_data> inis;
    std::vector<int> to_delete;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->RemoveData(recv_ptrn,to_delete);
    for(auto idx = to_delete.rbegin(); idx != to_delete.rend(); ++idx)
      destroyAnalysis(ans[*idx]);
    std::vector<int> to_add;
    std::vector<int> empty;
    size_t recv_init_ptrn = cs->AddData(recv_ptrn,empty,to_add);
    int ii = 0;
    ans.resize(ans.size()+to_add.size());
    cs->Communicate(recv_init_ptrn,
                    hdrs,
                    amsi::mpi_type<bio::micro_fo_header>());
    cs->Communicate(recv_init_ptrn,
                    prms,
                    amsi::mpi_type<bio::micro_fo_params>());
    cs->Communicate(recv_init_ptrn,
                    inis,
                    amsi::mpi_type<bio::micro_fo_init_data>());
    gmi_model * nl_mdl = gmi_load(".null");
    MPI_Comm slf;
    MPI_Comm_dup(MPI_COMM_SELF,&slf);
    PCU_Switch_Comm(slf);
    for(auto rve = ans.begin(); rve != ans.end(); ++rve)
    {
      if(*rve == NULL)
      {
        micro_fo_header & hdr = hdrs[ii];
        micro_fo_params & prm = prms[ii];
        micro_fo_init_data & dat = inis[ii];
        int tp = hdr.data[RVE_TYPE];
        int rnd = rand() % rve_tp_cnt[tp];
        apf::Mesh * msh_cpy = apf::createMdsMesh(nl_mdl,fns[tp][rnd]->msh);
        std::string fbr_rct_str("fiber_reaction");
        // copy fiber_reaction tag from origin mesh and set the reactions vector inside the fn
        amsi::copyIntTag(fbr_rct_str,fns[tp][rnd]->msh,msh_cpy,1,1);
        // copy ids for edges to make debugging easier
        AMSI_DEBUG(std::string id_tg_str("id"));
        AMSI_DEBUG(amsi::copyIntTag(id_tg_str,fns[tp][rnd]->msh,msh_cpy,1,1));
        FiberNetwork * fn = new FiberNetwork(msh_cpy);
        fn->getFiberReactions() = fns[tp][rnd]->rctns; // hate this, fix
        *rve = initFromMultiscale(fn,sprs[tp][rnd],bfrs,hdr,prm,dat);
        ii++;
      }
    }
    PCU_Switch_Comm(AMSI_COMM_SCALE);
    cs->CommPattern_UpdateInverted(recv_ptrn,send_ptrn);
    cs->CommPattern_Assemble(send_ptrn);
    cs->CommPattern_Reconcile(send_ptrn);
  }
  void MultiscaleRVEAnalysis::run()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    int sim_complete = 0;
    while(!sim_complete)
    {
      int step_complete = 0;
      while(!step_complete)
      {
        // migration
        if(macro_iter == 0)
          updateCoupling();
        std::vector<micro_fo_data> data;
        cs->Communicate(recv_ptrn,data,amsi::mpi_type<micro_fo_data>());
        std::vector<micro_fo_result> results(data.size());
        int ii = 0;
        for(auto rve = ans.begin(); rve != ans.end(); ++rve)
        {
          applyMultiscaleCoupling(*rve,&data[ii]);
          FiberRVEIteration itr(*rve);
          double prv_nrm = 1.0;
          auto val_gen = [&]() -> double
            {
              double nrm = (*rve)->ops->norm((*rve)->f);
              double val = fabs(prv_nrm - nrm);
              prv_nrm = nrm;
              return val;
            };
          auto eps_gen = [ ](int) -> double { return 1e-8; };
          auto ref_gen = [&]() -> double
            {
              return 1.0;
              /*
              static double nrm_f0 = 0.0;
              if(itr.iteration() == 1)
                nrm_f0 = (*rve)->ops->norm((*rve)->f0);
              return nrm_f0;
              */
            };
          amsi::UpdatingConvergence<decltype(&val_gen), decltype(&eps_gen), decltype(&ref_gen)> resid_cnvrg(&itr,&val_gen,&eps_gen,&ref_gen);
          amsi::Convergence * ptr[] = {&resid_cnvrg};
          amsi::MultiConvergence cnvrg(&ptr[0],&ptr[0]+1);
          // not a huge fan of this way vs adding it at the end of the multiconvergence, though this necessitates that the iteration reset happes at the END
          amsi::ResetIteration rst_iter(&cnvrg,&itr);
          amsi::numericalSolve(&itr,&cnvrg);
          // we've converged and have not reset the state of the vectors, matrices, and buffers
          // the inversion of the tangent stiffness matrix should be available in the buffers?
          recoverMultiscaleResults(*rve,&results[ii]);
          ++ii;
        }
        cs->Communicate(send_ptrn,results,amsi::mpi_type<micro_fo_result>());
        macro_iter++;
        cs->scaleBroadcast(M2m_id,&step_complete);
      }
      macro_iter = 0;
      macro_step++;
    }
  }
}
