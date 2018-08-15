#include "bioMultiscaleRVEAnalysis.h"
#include <amsiDetectOscillation.h>
#include <amsiNonlinearAnalysis.h>
#include <apfFEA.h>  // amsi
#include <apfMDS.h>
#include <apfMatrixUtil.h>    //amsi
#include <apfMeshIterator.h>  // amsi
#include <bioVerbosity.h>
#include <gmi.h>
#include <lasCSRCore.h>
#include <cassert>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOMultiscale.h"
#include "bioRVEVolumeTerms.h"
namespace bio
{
  void applyMultiscaleCoupling(FiberRVEAnalysis * ans, micro_fo_data * data)
  {
    int d = ans->rve->getDim();
    assert(d == 3 || d == 2);
    apf::Matrix3x3 F;
    for (int ei = 0; ei < d; ++ei)
      for (int ej = 0; ej < d; ++ej)
        F[ei][ej] = data->data[ei * d + ej];
    apf::Mesh * rve_msh = ans->rve->getMesh();
    apf::Field * rve_du = ans->rve->getdUField();
    apf::Field * rve_u = ans->rve->getUField();
    // apply deformation gradient to the rve cube
    ApplyDeformationGradient(F, rve_msh, rve_du, rve_u).run();
    apf::Mesh * fn_msh = ans->getFn()->getNetworkMesh();
    apf::Field * fn_du = ans->getFn()->getdUField();
    apf::Field * fn_u = ans->getFn()->getUField();
    // only apply first order continuation if we have derivative information
    // from the last iteration
    if (ans->dx_fn_dx_rve_set)
    {
      int dim = ans->rve->getDim();
      int nnd = ans->rve->numNodes();
      apf::DynamicVector rve_duv(nnd * dim);
      // array of the displacement field change since last itr of the rve cube
      apf::NewArray<apf::Vector3> rve_dus;
      apf::Element * du_elmt =
          apf::createElement(rve_du, ans->rve->getMeshEnt());
      apf::getVectorNodes(du_elmt, rve_dus);
      apf::destroyElement(du_elmt);
      apf::NewArray<int> rve_dofs;
      apf::getElementNumbers(
          ans->rve->getNumbering(), ans->rve->getMeshEnt(), rve_dofs);
      // reorder rve dofs with dofids since dx_fn_dx_rve uses them
      for (int nd = 0; nd < nnd; ++nd)
      {
        for (int dm = 0; dm < dim; ++dm)
        {
          rve_duv[rve_dofs[nd * dim + dm]] = rve_dus[nd][dm];
        }
      }
      int fn_dofs = ans->getFn()->getDofCount();
      apf::DynamicVector du_fn(fn_dofs);
      apf::multiply(ans->dx_fn_dx_rve, rve_duv, du_fn);
      ApplyDeformationGradient app_F(F, fn_msh, fn_du, fn_u);
      // Apply the deformation gradient to the boundary nodes and do not add
      // extra term from first order continuation e.g. set du_fn to zero on
      // boundary
      for (auto vrt = ans->bnd_nds[RVE::side::all].begin();
           vrt != ans->bnd_nds[RVE::side::all].end();
           ++vrt)
      {
        for (int dm = 0; dm < dim; ++dm)
        {
          int dof = apf::getNumber(ans->getFn()->getUNumbering(), *vrt, 0, dm);
          du_fn(dof) = 0.0;
        }
        app_F.inEntity(*vrt);
        app_F.atNode(0);
        app_F.outEntity();
      }
      // note that we want to accumulate the du*fx_fn_dx_rve onto the u field
      ApplySolution(fn_u, ans->getFn()->getUNumbering(), &du_fn[0], 0, true)
          .run();
    }
    else
    {
      // apply deformation gradient to the fiber network mesh
      ApplyDeformationGradient(F, fn_msh, fn_du, fn_u).run();
    }
  }
  void recoverMicroscaleStress(FiberRVEAnalysis * ans, double * stress)
  {
    auto ops = las::getLASOps<las::sparskit>();
    int dim = ans->getFn()->getNetworkMesh()->getDimension();
    apf::Field * xpu = ans->getFn()->getXpUField();
    apf::Vector3 xpu_val;
    apf::Vector3 f_val;
    apf::Numbering * num = ans->getFn()->getUNumbering();
    int dofs[3];
    double * f = NULL;
    ops->get(ans->getF(), f);
    apf::Matrix3x3 strs(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for (auto nd = ans->bnd_nds[RVE::all].begin();
         nd != ans->bnd_nds[RVE::all].end();
         ++nd)
    {
      apf::getVector(xpu, *nd, 0, xpu_val);
      for (int ii = 0; ii < ans->getFn()->getDim(); ++ii)
        dofs[ii] = apf::getNumber(num, *nd, 0, ii);
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
    amsi::mat2VoigtVec(dim, sym_strs, &stress[0]);
  }
  void convertStress(FiberRVEAnalysis * ans, double * stress)
  {
    int dim = ans->getFn()->getNetworkMesh()->getDimension();
    int sigma_length = dim == 3 ? 6 : 3;
    // convert to a macro-scale term
    double vol = ans->rve->measureDu();
    double scale_conversion = ans->multi->getScaleConversion();
    for (int ii = 0; ii < sigma_length; ++ii)
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
    dS_dx_rve.setSize(sigma_length, rve_dof_cnt);
    apf::DynamicVector col(sigma_length);
    for (int rve_dof = 0; rve_dof < rve_dof_cnt; ++rve_dof)
    {
      apf::DynamicVector sig(strs);
      ds_dx_rve.getColumn(rve_dof, col);
      col /= vol;
      sig *= (dV_dx_rve[rve_dof]);
      sig *= 1 / (vol * vol);
      col -= sig;
      col *= cnv;
      dS_dx_rve.setColumn(rve_dof, col);
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
    fn->getNetworkMesh()->getPoint(vrt, 0, crd);
    int vrt_cnt = dim == 3 ? 4 : 2;
    apf::MeshEntity * sd_vrts[vrt_cnt];
    rve->getMesh()->getDownward(sd_ent, 0, &sd_vrts[0]);
    apf::NewArray<int> rve_dof_ids(vrt_cnt * dim);
    apf::getElementNumbers(rve->getNumbering(), sd_ent, rve_dof_ids);
    apf::NewArray<int> fn_dof_ids(dim);
    apf::getElementNumbers(fn->getUNumbering(), vrt, fn_dof_ids);
    // iterate over the line/face (dep on dim)
    // assumes dim = 3
    // get the vert opposite the current vert on the face, which is what the
    //  value given by the calculation below is applicable to,
    //  we are basically calculating barycentric coordinates using
    //  the fact the RVE is unit cube as a huge shortcut rather than
    //  inverting the fe mapping
    int opp_vrt[] = {2, 3, 0, 1};
    for (int ci = 0; ci < vrt_cnt; ++ci)
    {
      apf::Vector3 ci_crd;
      rve->getMesh()->getPoint(sd_vrts[ci], 0, ci_crd);
      // this method of calculating the area/length associated
      //  with the vertex is only valid on the reference domain
      //  as it depends on the verts of a face/edge being both
      //  coplanar and axis-aligned
      apf::Vector3 spn = ci_crd - crd;
      double ai = 1.0;
      for (int ii = 0; ii < dim; ++ii)
        ai *= ii == plnr_dim ? 1.0 : spn[ii];
      ai = fabs(ai / a);
      // iterate over the dim
      for (int ej = 0; ej < dim; ++ej)
        dRdx_rve(fn_dof_ids[ej], rve_dof_ids[opp_vrt[ci] * dim + ej]) = ai;
    }
  }
  // wierd area term matrix
  // this is calculated on the reference domain and makes assumptions based on
  // that fact
  void calcdR_dx_rve(apf::DynamicMatrix & dRdx_rve, FiberRVEAnalysis * ans)
  {
    int dim = ans->rve->getDim();
    int fn_dof_cnt = ans->getFn()->getDofCount();
    int rve_dof_cnt = ans->rve->numNodes() * dim;
    dRdx_rve.setSize(fn_dof_cnt, rve_dof_cnt);
    dRdx_rve.zero();
    // this assumes 3d, use dim to change the sides we iterate over
    int sds[] = {RVE::side::top,
                 RVE::side::bot,
                 RVE::side::lft,
                 RVE::side::rgt,
                 RVE::side::bck,
                 RVE::side::frt};
    for (int sd_idx = 0; sd_idx != RVE::side::all; ++sd_idx)
    {
      int sd = sds[sd_idx];
      RVE::side rve_sd = static_cast<RVE::side>(sd);
      apf::MeshEntity * sd_ent = ans->rve->getSide(rve_sd);
      double a = apf::measure(ans->rve->getMesh(), sd_ent);
      // also assuming 3d
      int plnr_dim =
          (sd == RVE::side::rgt || sd == RVE::side::lft)
              ? 0
              : (sd == RVE::side::bot || sd == RVE::side::top) ? 1 : 2;
      for (auto vrt = ans->bnd_nds[sd].begin(); vrt != ans->bnd_nds[sd].end();
           ++vrt)
        calcdRdx_rve_term(
            dRdx_rve, ans->getFn(), ans->rve, sd_ent, *vrt, plnr_dim, a);
    }
  }
  // ans->k needs to be modified during calculation ofr dS_dx_fn prior to
  // calling this
  void calcdx_fn_dx_rve(apf::DynamicMatrix & dx_fn_dx_rve,
                        FiberRVEAnalysis * ans,
                        apf::DynamicMatrix & dR_dx_rve)
  {
    auto ops = las::getLASOps<las::sparskit>();
    int dim = ans->rve->getDim();
    int fn_dof_cnt = ans->getFn()->getDofCount();
    int rve_dof_cnt = ans->rve->numNodes() * dim;
    apf::DynamicVector f(fn_dof_cnt);
    apf::DynamicVector u(fn_dof_cnt);
    auto vb = las::getVecBuilder<las::sparskit>(0);
    las::Vec * skt_f = vb->create(fn_dof_cnt, LAS_IGNORE, MPI_COMM_SELF);
    las::Vec * skt_u = vb->create(fn_dof_cnt, LAS_IGNORE, MPI_COMM_SELF);
    las::Solve * iluslv = las::createSparskitLUSolve(ans->getSlv(), 0.0);
    las::Solve * qslv = las::createSparskitQuickLUSolve(iluslv, 0.0);
    dx_fn_dx_rve.setSize(fn_dof_cnt, rve_dof_cnt);
    double * fptr = NULL;
    double * uptr = NULL;
    for (int ii = 0; ii < rve_dof_cnt; ++ii)
    {
      las::Solve * slv = ii == 0 ? iluslv : qslv;
      dR_dx_rve.getColumn(ii, f);
      ops->get(skt_f, fptr);
      std::copy(f.begin(), f.end(), fptr);
      slv->solve(ans->getK(), skt_u, skt_f);
      ops->get(skt_u, uptr);
      std::copy(uptr, uptr + fn_dof_cnt, u.begin());
      dx_fn_dx_rve.setColumn(ii, u);
    }
    vb->destroy(skt_f);
    vb->destroy(skt_u);
    delete qslv;
    delete iluslv;
  }
  // setup additional solves with the rows of dR_dx_rve as the force vector
  // need to modify matrix as done in calc_precond in old code prior to these
  // solves, this also produces dS_dx_fn the change of the stresses on the
  // boundary of the RVE w.r.t. the fiber network coordinates
  void calcds_dx_fn(apf::DynamicMatrix & ds_dx_fn, FiberRVEAnalysis * ans)
  {
    auto ops = las::getLASOps<las::sparskit>();
    int dim = ans->rve->getDim();
    int sigma_length = dim == 3 ? 6 : 3;
    int fn_dof_cnt = ans->getFn()->getDofCount();
    ds_dx_fn.setSize(sigma_length, fn_dof_cnt);
    ds_dx_fn.zero();
    apf::Numbering * dofs = ans->getFn()->getUNumbering();
    double * F = NULL;
    ops->get(ans->getF(), F);
    // iterate over all fibers with a node on the boundary
    for (auto vrt = ans->bnd_nds[RVE::all].begin();
         vrt != ans->bnd_nds[RVE::all].end();
         ++vrt)
    {
      apf::Mesh * fn_msh = ans->getFn()->getNetworkMesh();
      assert(fn_msh->countUpward(*vrt) == 1);
      apf::Vector3 Nx;
      // need the displaced node coordinates, not the reference coords
      apf::getVector(ans->getFn()->getXpUField(), *vrt, 0, Nx);
      // fn_msh->getPoint(*vrt,0,Nx);
      apf::MeshEntity * edg = fn_msh->getUpward(*vrt, 0);
      apf::MeshEntity * vrts[2];
      fn_msh->getDownward(edg, 0, &vrts[0]);
      int bnd_dofs[dim];
      for (int dd = 0; dd < dim; ++dd)
        bnd_dofs[dd] = apf::getNumber(dofs, *vrt, 0, dd);
      for (int vrt = 0; vrt < 2; ++vrt)
      {
        for (int dd = 0; dd < dim; ++dd)
        {
          int dof = apf::getNumber(dofs, vrts[vrt], 0, dd);
          apf::Vector3 ks;
          ks.zero();
          for (int ii = 0; ii < dim; ++ii)
          {
            ks[ii] = las::getSparskitMatValue(ans->getK(), bnd_dofs[ii], dof);
            ks[ii] *= -1.0;
          }
          apf::Vector3 delta;
          for (int ii = 0; ii < dim; ++ii)
            delta[ii] = (dof == bnd_dofs[ii] ? 1.0 : 0.0);
          for (int i = dim; i < 3; ++i)
            delta[i] = -1;
          apf::Vector3 fs;
          for (int ii = 0; ii < dim; ++ii)
            fs[ii] = F[bnd_dofs[ii]];
          apf::Matrix3x3 m =
              apf::tensorProduct(ks, Nx) + apf::tensorProduct(delta, fs);
          apf::Matrix3x3 ms = amsi::symmetricPart(m);
          apf::DynamicVector vls(sigma_length);
          amsi::mat2VoigtVec(dim, ms, &vls(0));
          for (int ii = 0; ii < sigma_length; ++ii)
            ds_dx_fn(ii, dof) += vls(ii);
          for (int ii = 0; ii < dim; ++ii)
            las::setSparskitMatValue(ans->getK(), bnd_dofs[ii], dof, 0.0);
        }
      }
      for (int ii = 0; ii < dim; ++ii)
        las::setSparskitMatValue(ans->getK(), bnd_dofs[ii], bnd_dofs[ii], 1.0);
    }
  }
  void recoverStressDerivs(FiberRVEAnalysis * ans,
                           double * sigma,
                           double * dstrss_drve)
  {
    int dim = ans->rve->getDim();
    // int fn_dof_cnt = ans->getFn()->getDofCount();
    int sigma_length = dim == 3 ? 6 : 3;
    apf::DynamicMatrix dR_dx_rve;
    calcdR_dx_rve(dR_dx_rve, ans);
    apf::DynamicMatrix ds_dx_fn;
    calcds_dx_fn(ds_dx_fn, ans);
    // zero rows in the matrix w/ boundary conditions
    // this effects the force vector which is used in the calculation of
    // dS_dx_fn above so we must do it after.
    applyRVEBC(ans->bnd_nds[RVE::all].begin(),
               ans->bnd_nds[RVE::all].end(),
               ans->getFn()->getUNumbering(),
               ans->getK(),
               ans->getF());
    calcdx_fn_dx_rve(ans->dx_fn_dx_rve, ans, dR_dx_rve);
    apf::DynamicMatrix ds_dx_rve;
    apf::multiply(ds_dx_fn, ans->dx_fn_dx_rve, ds_dx_rve);
    ans->dx_fn_dx_rve_set = true;
    // calculate volume derivative
    apf::DynamicVector dV_dx_rve;
    CalcdV_dx_rve calcdv_dx_rve(
        2, ans->rve->getUField(), ans->rve->getNumbering());
    apf::MeshElement * mlm =
        apf::createMeshElement(ans->rve->getXpUField(), ans->rve->getMeshEnt());
    calcdv_dx_rve.process(mlm);
    calcdv_dx_rve.getdVdxrve(dV_dx_rve);
    apf::destroyMeshElement(mlm);
    apf::DynamicVector stress(sigma_length);
    memcpy(&stress[0], sigma, sizeof(double) * sigma_length);
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
    apf::multiply(dS_dx_rve, dx_rve_dx_fe, dS_dx);
    amsi::mat2Array(dS_dx, dstrss_drve);
  }
  void recoverMultiscaleResults(FiberRVEAnalysis * ans, micro_fo_result * data)
  {
    auto ops = las::getLASOps<las::sparskit>();
    // rebuild everything since we want the force vector without
    // boundary conditions anyway
    ops->zero(ans->getK());
    ops->zero(ans->getF());
    apf::Mesh * fn = ans->getFn()->getNetworkMesh();
    apf::MeshEntity * me = NULL;
    apf::MeshIterator * it = fn->begin(1);
    while ((me = fn->iterate(it)))
    {
      apf::MeshElement * mlm = apf::createMeshElement(fn, me);
      ans->es->process(mlm);
      apf::destroyMeshElement(mlm);
    }
    fn->end(it);
    double * S = &data->data[0];
    recoverMicroscaleStress(ans, S);
    double * Q = &data->data[6];
    recoverAvgVolStress(ans, Q);
    double * dS_dx_fe = &data->data[9];
    // calculate macroscale stress derivs
    // requires microscale stress
    recoverStressDerivs(ans, S, dS_dx_fe);
    // convert microscale stress to macroscale
    convertStress(ans, S);
  }
  void recoverMultiscaleStepResults(FiberRVEAnalysis * ans,
                                    micro_fo_header & hdr,
                                    micro_fo_params & prm,
                                    micro_fo_step_result * data)
  {
    double * ornt_3d = &data->data[0];
    double * ornt_2d = &data->data[9];
    double n[3];
    n[0] = prm.data[ORIENTATION_AXIS_X];
    n[1] = prm.data[ORIENTATION_AXIS_Y];
    n[2] = prm.data[ORIENTATION_AXIS_Z];
    if (hdr.data[COMPUTE_ORIENTATION_3D])
    {
      get3DOrientationTensor(ans->getFn(), ornt_3d);
    }
    if (hdr.data[COMPUTE_ORIENTATION_2D])
    {
      get2DOrientationTensor(ans->getFn(), n, ornt_2d);
    }
  }
  MultiscaleRVEAnalysis::~MultiscaleRVEAnalysis()
  {
    delete bfrs;
    for (auto v = vecs.begin(); v != vecs.end(); ++v)
    {
      delete (*v);
    }
    for (auto rve = ans.begin(); rve != ans.end(); ++rve)
    {
      destroyAnalysis(*rve);
      (*rve) = NULL;
    }
    assert(fns.size() == sprs.size() && fns.size() == dofs_cnt.size());
    // need to delete any fns that were not used as part of an analysis
    for (std::size_t i = 0; i < fns.size(); ++i)
    {
      for (int j = 0; j < rve_tp_cnt[i]; ++j)
      {
        delete fns[i][j];
        las::destroySparsity<las::CSR *>(sprs[i][j]);
        if (meshes[i][j])
        {
          meshes[i][j]->destroyNative();
          apf::destroyMesh(meshes[i][j]);
          meshes[i][j] = NULL;
        }
      }
      delete[] fns[i];
      delete[] sprs[i];
      delete[] dofs_cnt[i];
      delete[] meshes[i];
    }
    fns.clear();
    sprs.clear();
    dofs_cnt.clear();
    vecs.clear();
    meshes.clear();
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
      , hdrs()
      , prms()
      , slvr_prms()
      , slvr_int_prms()
      , sprs()
      , bfrs(NULL)
      , vecs()
      , dofs_cnt()
      , macro_iter(0)
      , macro_step(0)
  {
    M2m_id = amsi::getRelationID(amsi::getMultiscaleManager(),
                                 amsi::getScaleManager(),
                                 "macro",
                                 "micro_fo");
    m2M_id = amsi::getRelationID(amsi::getMultiscaleManager(),
                                 amsi::getScaleManager(),
                                 "micro_fo",
                                 "macro");
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
    rve_tp_lg = amsi::activateLog("rve_type_log");
  }
  void MultiscaleRVEAnalysis::initCoupling()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    rve_dd = amsi::createDataDistribution(amsi::getLocal(), "macro_fo_data");
    recv_ptrn = cs->RecvCommPattern("micro_fo_data", "macro", "macro_fo_data",
                                    "micro_fo");
    cs->CommPattern_Reconcile(recv_ptrn);
    send_ptrn = cs->CommPattern_UseInverted(recv_ptrn, "macro_fo_data",
                                            "micro_fo", "macro");
    cs->CommPattern_Assemble(send_ptrn);
    cs->CommPattern_Reconcile(send_ptrn);
  }
  void MultiscaleRVEAnalysis::initAnalysis()
  {
    int num_rve_tps = 0;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->scaleBroadcast(M2m_id, &num_rve_tps);
    char ** rve_tp_dirs = new char *[num_rve_tps];
    MPI_Request rqsts[num_rve_tps];
    // the order of receipt might be non-deterministic. need to handle that
    for (int ii = 0; ii < num_rve_tps; ++ii)
    {
      int cnt = 0;
      while ((cnt = cs->aRecvBroadcastSize<char>(M2m_id)) == 0)
      {
      }
      rve_tp_dirs[ii] = new char[cnt];
      // don't have to block to wait since we know the message was available for
      // size info
      cs->aRecvBroadcast(&rqsts[ii], M2m_id, &rve_tp_dirs[ii][0], cnt);
    }
    MPI_Status stss[num_rve_tps];
    MPI_Waitall(num_rve_tps, &rqsts[0], &stss[0]);
    MPI_Request hdr_rqst;
    rve_tp_cnt.resize(num_rve_tps);
    cs->aRecvBroadcast(&hdr_rqst, M2m_id, &rve_tp_cnt[0], num_rve_tps);
    MPI_Status hdr_sts;
    MPI_Waitall(1, &hdr_rqst, &hdr_sts);
    // Read in all the fiber network meshes and reactions
    for (int ii = 0; ii < num_rve_tps; ++ii)
    {
      fns.push_back(new FiberNetworkReactions *[rve_tp_cnt[ii]]);
      sprs.push_back(new las::Sparsity *[rve_tp_cnt[ii]]);
      dofs_cnt.push_back(new int[rve_tp_cnt[ii]]);
      meshes.push_back(new apf::Mesh2 *[rve_tp_cnt[ii]]);
    }
    int dof_max = -1;
    PCU_Switch_Comm(MPI_COMM_SELF);
    for (int ii = 0; ii < num_rve_tps; ii++)
    {
      for (int jj = 0; jj < rve_tp_cnt[ii]; ++jj)
      {
        std::stringstream fl;
        fl << rve_tp_dirs[ii] << jj + 1 << ".txt";
        FiberNetworkReactions * fn_rctns = new FiberNetworkReactions;
        meshes[ii][jj] = loadFromFile(fl.str());
        apf::Mesh2 * fn = meshes[ii][jj];
        fn_rctns->fileName = fl.str();
        fl << ".params";
        loadParamsFromFile(fn, fl.str(), std::back_inserter(fn_rctns->rctns));
        fns[ii][jj] = fn_rctns;
        apf::Field * u = apf::createLagrangeField(fn, "u", apf::VECTOR, 1);
        apf::zeroField(u);
        apf::Numbering * n = apf::createNumbering(u);
        int dofs = apf::NaiveOrder(n);
        sprs[ii][jj] = las::createCSR(n, dofs);
        dofs_cnt[ii][jj] = dofs;
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
    std::vector<micro_fo_init_data> inis;
    std::vector<int> to_delete;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->RemoveData(recv_ptrn, to_delete);
    for (auto idx = to_delete.rbegin(); idx != to_delete.rend(); ++idx)
    {
      // FIXME need to clear out LinearStructs,
      // sparsity, dofs_cnt, etc here.
      // We also cannot just destroy the analysis w/o removing the corresponding
      // terms from the ans vector see vector erase
      // destroyAnalysis(ans[*idx]);
    }
    std::vector<int> to_add;
    std::vector<int> empty;
    size_t recv_init_ptrn = cs->AddData(recv_ptrn, empty, to_add);
    ans.resize(ans.size() + to_add.size());
    cs->Communicate(recv_init_ptrn, hdrs,
                    amsi::mpi_type<bio::micro_fo_header>());
    cs->Communicate(recv_init_ptrn, prms,
                    amsi::mpi_type<bio::micro_fo_params>());
    cs->Communicate(recv_init_ptrn, inis,
                    amsi::mpi_type<bio::micro_fo_init_data>());
    cs->Communicate(recv_init_ptrn, slvr_prms,
                    amsi::mpi_type<bio::micro_fo_solver>());
    cs->Communicate(recv_init_ptrn, slvr_int_prms,
                    amsi::mpi_type<bio::micro_fo_int_solver>());
    PCU_Switch_Comm(MPI_COMM_SELF);
    int ii = 0;
    for (auto rve = ans.begin(); rve != ans.end(); ++rve)
    {
      if (*rve == NULL)
      {
        micro_fo_header & hdr = hdrs[ii];
        micro_fo_params & prm = prms[ii];
        micro_fo_init_data & dat = inis[ii];
        micro_fo_solver & slvr_prm = slvr_prms[ii];
        micro_fo_int_solver & slvr_int_prm = slvr_int_prms[ii];
        int tp = hdr.data[RVE_TYPE];
        int rnd = rand() % rve_tp_cnt[tp];
        apf::Mesh * msh_cpy =
            apf::createMdsMesh(gmi_load(".null"), meshes[tp][rnd]);
        FiberNetwork * fn = new FiberNetwork(msh_cpy);
        fn->setFiberReactions(fns[tp][rnd]->rctns);
        vecs.push_back(
            createLinearStructs(dofs_cnt[tp][rnd], sprs[tp][rnd], bfrs));
        *rve = initFromMultiscale(
            fn, vecs[ii], hdr, prm, dat, slvr_prm, slvr_int_prm);
        fn->setRVEType(ii);
        BIO_V2(
            // print the list of fiber network names to file
            amsi::log(rve_tp_lg)
                << ii << " " << fns[tp][rnd]->fileName << std::endl;)
        ++ii;
      }
    }
    PCU_Switch_Comm(AMSI_COMM_SCALE);
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::ofstream rve_tp_lg_fs(
        amsi::fs->getResultsDir() + "/rve_tp." + std::to_string(rank) + ".log",
        std::ios::out | std::ios::app);
    amsi::flushToStream(rve_tp_lg, rve_tp_lg_fs);
    cs->CommPattern_UpdateInverted(recv_ptrn, send_ptrn);
    cs->CommPattern_Assemble(send_ptrn);
    cs->CommPattern_Reconcile(send_ptrn);
  }
  struct val_gen
  {
    val_gen(FiberRVEAnalysis * a) : an(a), prv_nrm(1.0) {}
    FiberRVEAnalysis * an;
    double prv_nrm;
    double operator()()
    {
      auto ops = las::getLASOps<las::sparskit>();
      double nrm = ops->norm(an->getF());
      double val = fabs(prv_nrm - nrm);
      prv_nrm = nrm;
      return val;
    }
  };
  struct eps_gen
  {
    eps_gen(double eps) : eps(eps) {}
    double operator()(int) { return eps; }
    protected:
    double eps;
  };
  struct ref_gen
  {
    double operator()() { return 1.0; }
  };
  void MultiscaleRVEAnalysis::run()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    bool sim_complete = false;
    while (!sim_complete)
    {
      bool step_complete = false;
      while (!step_complete)
      {
        // migration
        if (macro_iter == 0) updateCoupling();
        std::vector<micro_fo_data> data;
        cs->Communicate(recv_ptrn, data, amsi::mpi_type<micro_fo_data>());
        std::vector<micro_fo_result> results(data.size());
        PCU_Switch_Comm(MPI_COMM_SELF);
        int rank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int ii = 0;
        BIO_V1(double t0 = MPI_Wtime();)
        FiberRVEAnalysis * tmpRVE = NULL;
        for (auto rve = ans.begin(); rve != ans.end(); ++rve)
        {
          unsigned int maxMicroAttempts = 0;   // parameter
          unsigned int microAttemptCutFactor;  // parameter
          bool solveSuccess = false;
          unsigned int microAttemptCount = 1;
          unsigned int attemptCutFactor;
          do
          {
            // create a deep copy of the analysis
            // Note the current implementation of copy does not deep copy
            // the sparskit matrices, vectors, or solver
            tmpRVE = copyAnalysis(*rve);
            val_gen vg(tmpRVE);
            eps_gen eg(tmpRVE->solver_eps);
            ref_gen rg;
            maxMicroAttempts = tmpRVE->max_cut_attempt;
            microAttemptCutFactor = tmpRVE->attempt_cut_factor;
            attemptCutFactor =
                std::pow(microAttemptCutFactor, microAttemptCount - 1);
            BIO_V1(if (attemptCutFactor > 1) std::cout
                       << "Micro Attempt: " << microAttemptCount
                       << " cutting the original deformation gradient by: "
                       << attemptCutFactor << " on rank: " << rank << "\n";)
            assert(maxMicroAttempts > 0);
            micro_fo_data appliedDefm;
            bool microIterSolveSuccess = true;
            for (unsigned int microAttemptIter = 1;
                 microAttemptIter <= attemptCutFactor;
                 ++microAttemptIter)
            {
              BIO_V3(std::cout << "Rank: " << rank << " F=";)
              for (int j = 0; j < 9; ++j)
              {
                appliedDefm.data[j] =
                    (data[ii].data[j] * microAttemptIter) / attemptCutFactor;
                BIO_V3(std::cout << appliedDefm.data[j] << " ";)
              }
              BIO_V3(std::cout << "\n";)
              applyMultiscaleCoupling(tmpRVE, &appliedDefm);
              FiberRVEIteration rveItr(tmpRVE);
              std::vector<amsi::Iteration *> itr_stps = {&rveItr};
              amsi::MultiIteration itr(itr_stps.begin(), itr_stps.end());
              amsi::UpdatingConvergence<decltype(&vg),
                                        decltype(&eg),
                                        decltype(&rg)>
                  resid_cnvrg(&itr, &vg, &eg, &rg);
              std::vector<amsi::Convergence *> cnvrg_stps = {&resid_cnvrg};
              amsi::MultiConvergence cnvrg(cnvrg_stps.begin(),
                                           cnvrg_stps.end());
              amsi::Iteration * osc_itr =
                  amsi::createOscillationDetection<decltype(&resid_cnvrg)>(
                      tmpRVE->detect_osc_type,
                      &resid_cnvrg,
                      &rveItr,
                      tmpRVE->max_itrs,
                      tmpRVE->prev_itr_factor);
              itr.addIteration(osc_itr);
              // solve is successful if the current solve and all previous
              // cutIterations are successful
              microIterSolveSuccess =
                  (amsi::numericalSolve(&itr, &cnvrg) && microIterSolveSuccess);
              // cleanup the oscillation detection memory
              delete osc_itr;
              // don't bother computing the rest of the attempt if any
              // subiteration fails, for our current use case we don't care what
              // made us fail, we will try to reduce the load and try again.
              if (!microIterSolveSuccess) break;
            }
            // if the attempt was completely successful then the overall solve
            // was successful
            if (microIterSolveSuccess)
            {
              solveSuccess = true;
              destroyAnalysis(*rve);
              (*rve) = tmpRVE;
              tmpRVE = NULL;
            }
            ++microAttemptCount;
          } while (solveSuccess == false &&
                   (microAttemptCount <= maxMicroAttempts));
          if (!solveSuccess)
          {
            std::cerr << "RVE: " << (*rve)->getFn()->getRVEType()
                      << " failed to converge in " << microAttemptCount - 1
                      << " attempts on processor " << rank << std::endl;
            std::abort();  // should I use MPI_Abort() here?
          }
          // we've converged and have not reset the state of the vectors,
          // matrices, and buffers the inversion of the tangent stiffness matrix
          // should be available in the buffers?
          recoverMultiscaleResults(*rve, &results[ii]);
          ii++;
#ifdef WRITE_MICRO_PER_ITER
          std::stringstream sout;
          int rnk = -1;
          MPI_Comm_rank(AMSI_COMM_SCALE, &rnk);
          sout << "rnk_" << rnk << "_fn_" << ii << "_step_" << macro_step
               << "_iter_" << macro_iter;
          apf::writeVtkFiles(
              sout.str().c_str(), (*rve)->fn->getNetworkMesh(), 1);
#endif
        }
        BIO_V1(double t1 = MPI_Wtime();)
        BIO_V1(std::cout << "Computed " << ans.size() << " RVEs in " << t0 - t1
                         << " seconds." << std::endl;)
        PCU_Switch_Comm(AMSI_COMM_SCALE);
        cs->Communicate(send_ptrn, results, amsi::mpi_type<micro_fo_result>());
        macro_iter++;
        cs->scaleBroadcast(M2m_id, &step_complete);
      }
      // get the size of the step results vector
      std::vector<micro_fo_step_result> step_results(hdrs.size());
      // recover step results and set the step results vector
      int i = 0;
      PCU_Switch_Comm(MPI_COMM_SELF);
      for (auto rve = ans.begin(); rve != ans.end(); ++rve)
      {
        micro_fo_header & hdr = hdrs[i];
        micro_fo_params & prm = prms[i];
        recoverMultiscaleStepResults(*rve, hdr, prm, &step_results[i]);
        ++i;
      }
      PCU_Switch_Comm(AMSI_COMM_SCALE);
      // communicate the step results back to the macro scale
      amsi::ControlService * cs = amsi::ControlService::Instance();
      cs->Communicate(
          send_ptrn, step_results, amsi::mpi_type<micro_fo_step_result>());
#ifdef WRITE_MICRO_PER_STEP
      for (auto rve = ans.begin(); rve != ans.end(); ++rve)
      {
        std::stringstream sout;
        int rnk = -1;
        MPI_Comm_rank(AMSI_COMM_SCALE, &rnk);
        int ii = 0;
        sout << "rnk_" << rnk << "_fn_" << ii << "_step_" << macro_step;
        apf::writeVtkFiles(sout.str().c_str(), (*rve)->fn->getNetworkMesh(), 1);
        sout.str("");
        ii++;
      }
#endif
      macro_iter = 0;
      macro_step++;
      cs->scaleBroadcast(M2m_id, &sim_complete);
    }
  }
}  // namespace bio
