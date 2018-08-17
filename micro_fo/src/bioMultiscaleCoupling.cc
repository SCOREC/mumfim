#include "bioMultiscaleCoupling.h"
#include <apf.h>
#include <apfMatrixUtil.h>
#include <cassert>
#include "bioFiberRVEAnalysis.h"
#include "bioMultiscaleMicroFOParams.h"
#include "bioRVEVolumeTerms.h"
#include "bioUtil.h"
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
}  // namespace bio
