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
    if(ans->dx_fn_dx_rve_set == true)
    {
      int dim = ans->rve->getDim();
      int nnd = ans->rve->numNodes();
      apf::DynamicVector rve_uv(nnd*dim);
      apf::NewArray<apf::Vector3> rve_us;
      apf::getVectorNodes(ans->rve->getElement(),rve_us);
      apf::NewArray<int> rve_dofs;
      apf::getElementNumbers(ans->rve->getNumbering(),ans->rve->getMeshEnt(),rve_dofs);
      // reorder the rve dofs using dofids because
      // dx_fn_dx_rve uses them
      for(int nd = 0; nd < nnd; ++nd)
        for(int dm = 0; dm < dim; ++dm)
          rve_uv[rve_dofs[nd*dim+dm]] = rve_us[nd][dm];
      int fn_dofs = ans->fn->getDofCount();
      apf::DynamicVector u_fn(fn_dofs);
      apf::multiply(ans->dx_fn_dx_rve,rve_uv,u_fn);
      ApplyDeformationGradient app_F(F,fn_msh,fn_du,fn_u);
      for(auto vrt = ans->bnd_nds[RVE::side::all].begin(); vrt != ans->bnd_nds[RVE::side::all].end(); ++vrt)
      {
        for(int dm = 0; dm < dim; ++dm)
        {
          int dof = apf::getNumber(ans->fn->getUNumbering(),
                                   *vrt,0,dm);
          u_fn(dof) = 0.0;
        }
        app_F.inEntity(*vrt);
        app_F.atNode(0);
        app_F.outEntity();
      }
      // only changes 'u', 'du' is changed for boundary nodes but not internal nodes, will this cause
      // any issues?
      //
      ApplySolution(ans->fn->getUNumbering(),&u_fn[0],0,true).apply(ans->fn->getUField());
      apf::zeroField(ans->fn->getdUField());
      /*
      if(!PCU_Comm_Self())
      {
        std::ofstream fout("u_fld");
        amsi::PrintField(ans->fn->getXpUField(),fout).run();
      }
      */
    }
    else
      ApplyDeformationGradient(F,fn_msh,fn_du,fn_u).run();
  }
  // might be duplicated from bioFiberRVEAnalysis.cc ?
  void recoverMicroscaleStress(FiberRVEAnalysis * ans, double * stress)
  {
    int dim = ans->fn->getNetworkMesh()->getDimension();
    apf::Field * xpu = ans->fn->getXpUField();
    apf::Vector3 xpu_val;
    apf::Vector3 f_val;
    apf::Numbering * num = ans->fn->getUNumbering();
    int dofs[3] {};
    double * f = NULL;
    ans->ops->get(ans->f,f);
    apf::Matrix3x3 strs(0.0,0.0,0.0,
                        0.0,0.0,0.0,
                        0.0,0.0,0.0);
    for(auto nd = ans->bnd_nds[RVE::all].begin(); nd != ans->bnd_nds[RVE::all].end(); ++nd)
    {
      apf::getVector(xpu,*nd,0,xpu_val);
      for(int ii = 0; ii < ans->fn->getDim(); ++ii)
        dofs[ii] = apf::getNumber(num,*nd,0,ii);
      strs[0][0] += xpu_val[0] * f[dofs[0]];
      strs[0][1] += xpu_val[1] * f[dofs[0]];
      strs[0][2] += xpu_val[2] * f[dofs[0]];
      strs[1][0] += xpu_val[0] * f[dofs[1]];
      strs[1][1] += xpu_val[1] * f[dofs[1]];
      strs[1][2] += xpu_val[2] * f[dofs[1]];
      strs[2][0] += xpu_val[0] * f[dofs[2]];
      strs[2][1] += xpu_val[1] * f[dofs[2]];
      strs[2][2] += xpu_val[2] * f[dofs[2]];
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
  void dsdxrve_2_dSdxrve(const apf::DynamicMatrix & ds_dx_rve,
                         const apf::DynamicVector & strs,
                         const apf::DynamicVector & dV_dx_rve,
                         double vol,
                         double cnv,
                         apf::DynamicMatrix & dS_dx_rve)
  {
    int rve_dof_cnt = ds_dx_rve.getColumns();
    int sigma_length = ds_dx_rve.getRows();
    dS_dx_rve.setSize(sigma_length,rve_dof_cnt);
    apf::DynamicVector col(sigma_length);
    for(int rve_dof = 0; rve_dof < rve_dof_cnt; ++rve_dof)
    {
      apf::DynamicVector sig(strs);
      ds_dx_rve.getColumn(rve_dof,col);
      col /= vol;
      sig *= (dV_dx_rve[rve_dof]);
      sig *= 1 / (vol * vol);
      col -= sig;
      col *= cnv;
      dS_dx_rve.setColumn(rve_dof,col);
    }
  }
  void calcdRdx_rve_term(apf::DynamicMatrix & dRdx_rve,
                         FiberNetwork * fn,
                         RVE * rve,
                         apf::MeshEntity * sd_ent,
                         apf::MeshEntity * vrt,
                         int plnr_dim,
                         double a)
  {
    int dim = rve->getDim();
    apf::Vector3 crd;
    fn->getNetworkMesh()->getPoint(vrt,0,crd);
    int vrt_cnt = dim == 3 ? 4 : 2;
    apf::MeshEntity * sd_vrts[vrt_cnt];
    rve->getMesh()->getDownward(sd_ent,0,&sd_vrts[0]);
    apf::NewArray<int> rve_dof_ids(vrt_cnt * dim);
    apf::getElementNumbers(rve->getNumbering(),sd_ent,rve_dof_ids);
    apf::NewArray<int> fn_dof_ids(dim);
    apf::getElementNumbers(fn->getUNumbering(),vrt,fn_dof_ids);
    // iterate over the line/face (dep on dim)
    // assumes dim = 3
    // get the vert opposite the current vert on the face, which is what the
    //  value given by the calculation below is applicable to,
    //  we are basically calculating barycentric coordinates using
    //  the fact the RVE is unit cube as a huge shortcut rather than
    //  inverting the fe mapping
    int opp_vrt[] = {2,3,0,1};
    for(int ci = 0; ci < vrt_cnt; ++ci)
    {
      apf::Vector3 ci_crd;
      rve->getMesh()->getPoint(sd_vrts[ci],0,ci_crd);
      // this method of calculating the area/length associated
      //  with the vertex is only valid on the reference domain
      //  as it depends on the verts of a face/edge being both
      //  coplanar and axis-aligned
      apf::Vector3 spn = ci_crd - crd;
      double ai = 1.0;
      for(int ii = 0; ii < dim; ++ii)
        ai *= ii == plnr_dim ? 1.0 : spn[ii];
      ai = fabs(ai / a);
      // iterate over the dim
      for(int ej = 0; ej < dim; ++ej)
        dRdx_rve(fn_dof_ids[ej],rve_dof_ids[opp_vrt[ci]*dim + ej]) = ai;
    }
  }
  // wierd area term matrix
  // this is calculated on the reference domain and makes assumptions based on that fact
  void calcdR_dx_rve(apf::DynamicMatrix & dRdx_rve,
                    FiberRVEAnalysis * ans)
  {
    int dim = ans->rve->getDim();
    int fn_dof_cnt = ans->fn->getDofCount();
    int rve_dof_cnt = ans->rve->numNodes()*dim;
    dRdx_rve.setSize(fn_dof_cnt,rve_dof_cnt);
    dRdx_rve.zero();
    // this assumes 3d, use dim to change the sides we iterate over
    int sds[] = {RVE::side::top, RVE::side::bot, RVE::side::lft, RVE::side::rgt, RVE::side::bck, RVE::side::frt};
    for(int sd_idx = 0; sd_idx != RVE::side::all; ++sd_idx)
    {
      int sd = sds[sd_idx];
      RVE::side rve_sd = static_cast<RVE::side>(sd);
      apf::MeshEntity * sd_ent = ans->rve->getSide(rve_sd);
      double a = apf::measure(ans->rve->getMesh(),sd_ent);
      // also assuming 3d
      int plnr_dim = (sd == RVE::side::rgt || sd == RVE::side::lft) ? 0 :
        (sd == RVE::side::bot || sd == RVE::side::top) ? 1 : 2;
      for(auto vrt = ans->bnd_nds[sd].begin(); vrt != ans->bnd_nds[sd].end(); ++vrt)
        calcdRdx_rve_term(dRdx_rve,ans->fn,ans->rve,sd_ent,*vrt,plnr_dim,a);
    }
  }
  // ans->k needs to be modified during calculation ofr dS_dx_fn prior to calling this
  void calcdx_fn_dx_rve(apf::DynamicMatrix & dx_fn_dx_rve,
                        FiberRVEAnalysis * ans,
                        apf::DynamicMatrix & dR_dx_rve)
  {
    int dim = ans->rve->getDim();
    int fn_dof_cnt = ans->fn->getDofCount();
    int rve_dof_cnt = ans->rve->numNodes()*dim;
    apf::DynamicVector f(fn_dof_cnt);
    apf::DynamicVector u(fn_dof_cnt);
    las::Vec * skt_f = las::createSparskitVector(fn_dof_cnt);
    las::Vec * skt_u = las::createSparskitVector(fn_dof_cnt);
    las::LasSolve * slv = las::createSparskitLUSolve(las::getSparskitBuffers(ans->slv),0.0);
    dx_fn_dx_rve.setSize(fn_dof_cnt,rve_dof_cnt);
    double * fptr = NULL;
    double * uptr = NULL;
    for(int ii = 0; ii < rve_dof_cnt; ++ii)
    {
      dR_dx_rve.getColumn(ii,f);
      ans->ops->get(skt_f,fptr);
      std::copy(f.begin(),f.end(),fptr);
      slv->solve(ans->k,skt_u,skt_f);
      ans->ops->get(skt_u,uptr);
      std::copy(uptr,uptr+fn_dof_cnt,u.begin());
      dx_fn_dx_rve.setColumn(ii,u);
    }
  }
  // setup additional solves with the rows of dR_dx_rve as the force vector
  // need to modify matrix as done in calc_precond in old code prior to these solves, this also produces dS_dx_fn the change of the stresses on the boundary of the RVE w.r.t. the fiber network coordinates
  void calcds_dx_fn(apf::DynamicMatrix & ds_dx_fn,
                    FiberRVEAnalysis * ans)
  {
    int dim = ans->rve->getDim();
    int sigma_length = dim == 3 ? 6 : 3;
    int fn_dof_cnt = ans->fn->getDofCount();
    ds_dx_fn.setSize(sigma_length,fn_dof_cnt);
    ds_dx_fn.zero();
    apf::Numbering * dofs = ans->fn->getUNumbering();
    double * F = NULL;
    ans->ops->get(ans->f,F);
    // iterate over all fibers with a node on the boundary
    for(auto vrt = ans->bnd_nds[RVE::all].begin(); vrt != ans->bnd_nds[RVE::all].end(); ++vrt)
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
            ds_dx_fn(ii,dof) += vls(ii);
          for(int ii = 0; ii < dim; ++ii)
            las::setSparskitMatValue(ans->k,bnd_dofs[ii],dof,0.0);
        }
      }
      for(int ii = 0; ii < dim; ++ii)
        las::setSparskitMatValue(ans->k,bnd_dofs[ii],bnd_dofs[ii],1.0);
    }
  }
  void recoverStressDerivs(FiberRVEAnalysis * ans, double * sigma, double * dstrss_drve)
  {
    int dim = ans->rve->getDim();
    //int fn_dof_cnt = ans->fn->getDofCount();
    int sigma_length = dim == 3 ? 6 : 3;
    apf::DynamicMatrix dR_dx_rve;
    calcdR_dx_rve(dR_dx_rve,ans);
    apf::DynamicMatrix ds_dx_fn;
    calcds_dx_fn(ds_dx_fn,ans);
    // zero rows in the matrix w/ boundary conditions
    // this effects the force vector which is used in the calculation of
    // dS_dx_fn above so we must do it after.
    applyRVEBC(ans->bnd_nds[RVE::all].begin(),ans->bnd_nds[RVE::all].end(),
               ans->fn->getUNumbering(),ans->ops,ans->k,ans->f);
    // calculate dx_fn_dx_rve;
    calcdx_fn_dx_rve(ans->dx_fn_dx_rve,ans,dR_dx_rve);
    apf::DynamicMatrix ds_dx_rve;
    apf::multiply(ds_dx_fn,ans->dx_fn_dx_rve,ds_dx_rve);
    ans->dx_fn_dx_rve_set = true;
    // calculate volume derivative
    apf::DynamicVector dV_dx_rve;
    CalcdV_dx_rve calcdv_dx_rve(2,
                                ans->rve->getUField(),
                                ans->rve->getNumbering());
    apf::MeshElement * mlm = apf::createMeshElement(ans->rve->getXpUField(),ans->rve->getMeshEnt());
    calcdv_dx_rve.process(mlm);
    calcdv_dx_rve.getdVdxrve(dV_dx_rve);
    apf::destroyMeshElement(mlm);
    apf::DynamicVector stress(sigma_length);
    memcpy(&stress[0],sigma,sizeof(double)*sigma_length);
    apf::DynamicMatrix dS_dx_rve;
    dsdxrve_2_dSdxrve(ds_dx_rve,
                      stress,
                      dV_dx_rve,
                      ans->rve->measureDu(),
                      ans->multi->getScaleConversion(),
                      dS_dx_rve);
    apf::DynamicMatrix dx_rve_dx_fe;
    ans->multi->calcdRVEdFE(dx_rve_dx_fe);
    apf::DynamicMatrix dS_dx;
    apf::multiply(dS_dx_rve,dx_rve_dx_fe,dS_dx);
    amsi::mat2Array(dS_dx,dstrss_drve);
  }
  void recoverMultiscaleResults(FiberRVEAnalysis * ans, micro_fo_result * data)
  {
    // rebuild everything since we want the force vector without
    // boundary conditions anyway
    ans->ops->zero(ans->k);
    ans->ops->zero(ans->f);
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
    PCU_Switch_Comm(MPI_COMM_SELF);
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
    PCU_Switch_Comm(AMSI_COMM_SCALE);
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
    PCU_Switch_Comm(MPI_COMM_SELF);
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
        PCU_Switch_Comm(MPI_COMM_SELF);
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
          auto eps_gen = [ ](int) -> double { return 1e-6; };
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
          ii++;
        }
        PCU_Switch_Comm(AMSI_COMM_SCALE);
        cs->Communicate(send_ptrn,results,amsi::mpi_type<micro_fo_result>());
        macro_iter++;
        cs->scaleBroadcast(M2m_id,&step_complete);
      }
      macro_iter = 0;
      macro_step++;
    }
  }
}
