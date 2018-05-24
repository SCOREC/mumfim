#include "bioRVE.h"
#include "bioFiberNetwork.h"
#include "lasSparskit.h"
#include <apfFunctions.h>
#include <apfMeshIterator.h>
#include <apfMeshUtil.h>
namespace bio
{
  template <typename O>
  void originCenterCube(O cbe_crnrs, double crd)
  {
    // APF orders differently than canonical
    // iosparametric hexs
    // front
    *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd * -1.0,crd * -1.0);
    *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd * -1.0,crd * -1.0);
    *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd *  1.0,crd * -1.0);
    *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd *  1.0,crd * -1.0);
    // back
    *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd * -1.0,crd *  1.0);
    *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd * -1.0,crd *  1.0);
    *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd *  1.0,crd *  1.0);
    *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd *  1.0,crd *  1.0);
  }
  template <typename O>
  void originCenterSquare(O cbe_crnrs, double crd)
  {
    // bottom
    *cbe_crnrs++ = apf::Vector3(-crd,-crd, 0.0);
    *cbe_crnrs++ = apf::Vector3( crd,-crd, 0.0);
    // top
    *cbe_crnrs++ = apf::Vector3( crd, crd, 0.0);
    *cbe_crnrs++ = apf::Vector3(-crd, crd, 0.0);
  }
  RVE::RVE(double cr, int d)
  : dim(d)
  , crd(cr)
  , cbe(NULL)
  , cbe_e(NULL)
  , cbe_u_e(NULL)
  , cbe_u(NULL)
  , cbe_du(NULL)
  , cbe_dof(NULL)
  {
    assert(d == 2 || d == 3);
    std::vector<apf::Vector3> cbe_crnrs;
    apf::Mesh::Type tp = apf::Mesh::TYPES;
    if(dim == 3)
    {
      originCenterCube(std::back_inserter(cbe_crnrs),crd);
      tp = apf::Mesh::HEX;
    }
    else if(dim == 2)
    {
      originCenterSquare(std::back_inserter(cbe_crnrs),crd);
      tp = apf::Mesh::QUAD;
    }
    cbe = amsi::makeSingleEntityMesh(tp,&cbe_crnrs[0]);
    apf::MeshIterator * it = cbe->begin(3);
    cbe_e = cbe->iterate(it);
    cbe->end(it);
    cbe_u = apf::createLagrangeField(cbe,"u",apf::VECTOR,1);
    apf::zeroField(cbe_u);
    cbe_du = apf::createLagrangeField(cbe,"du",apf::VECTOR,1);
    apf::zeroField(cbe_du);
    cbe_dof = apf::createNumbering(cbe_u);
    xpufnc = new amsi::XpYFunc(cbe->getCoordinateField(),cbe_u);
    cbe_xpu = apf::createUserField(cbe,"xpu",apf::VECTOR,apf::getShape(cbe_u),xpufnc);
    int cbe_dof_cnt = apf::NaiveOrder(cbe_dof);
    assert(dim == 3 ? cbe_dof_cnt == 24 : cbe_dof_cnt == 8);
    cbe_u_e = apf::createElement(cbe_u,cbe_e);
  }
  RVE::~RVE()
  {
    delete xpufnc;
  }
  apf::MeshEntity * RVE::getSide(side sd)
  {
    apf::MeshEntity * rslt = NULL;
    int plnr_dim =
      (sd == rgt || sd == lft) ? 0 :
      (sd == bot || sd == top) ? 1 : 2;
    double plnr_crd = sideCoord(sd);
    for(auto fc_itr = amsi::apfMeshIterator(cbe,2); fc_itr != amsi::apfEndIterator(cbe); ++fc_itr)
    {
      apf::MeshEntity * vrts[4];
      cbe->getDownward(*fc_itr,0,&vrts[0]);
      bool found = true;
      for(int ii = 0; ii < 4 && found; ++ii)
      {
        apf::Vector3 crd;
        cbe->getPoint(vrts[ii],0,crd);
        if(crd[plnr_dim] != plnr_crd)
          found = false;
      }
      if(found)
      {
        rslt = *fc_itr;
        break;
      }
    }
    return rslt;
  }
  void calcGlobalRVECoords(apf::DynamicArray<apf::Vector3> & rve_crds,
                           double rve_dim,
                           const apf::Vector3 & gbl_gss)
  {
    static double op[8][3] = {{-1.0,-1.0, 1.0},
                              { 1.0,-1.0, 1.0},
                              {-1.0,-1.0,-1.0},
                              { 1.0 -1.0,-1.0},
                              {-1.0, 1.0, 1.0},
                              { 1.0, 1.0, 1.0},
                              {-1.0, 1.0,-1.0},
                              { 1.0, 1.0,-1.0}};
    int sz = rve_crds.getSize();
    int d = sz == 8 ? 3 : 2;
    int o = d == 3 ? 0 : 2;
    for(int ii = 0; ii < sz; ii++)
      for(int jj = 0; jj < d; jj++)
        rve_crds[ii][jj] = gbl_gss[jj] + rve_dim*op[ii+o][jj];
  }
  void displaceRVE(RVE * rve,const apf::DynamicVector & du)
  {
    ApplySolution(rve->getNumbering(),&du[0],0,true).apply(rve->getUField());
  }
  void getRVEDisplacement(RVE * rve, apf::DynamicVector & u)
  {
    int ndofs = rve->getDim() * rve->numNodes();
    int sz = u.getSize();
    if(sz != ndofs)
      u.setSize(ndofs);
    amsi::WriteOp wrt;
    amsi::ToArray(rve->getNumbering(),rve->getUField(),&u(0),0,&wrt).run();
  }
  void getRVEReferenceCoords(RVE * rve, apf::DynamicVector & xyz_0)
  {
    int ndofs = rve->getDim() * rve->numNodes();
    int sz = xyz_0.getSize();
    if(sz != ndofs)
      xyz_0.setSize(ndofs);
    amsi::WriteOp wrt;
    amsi::ToArray(rve->getNumbering(),rve->getMesh()->getCoordinateField(),&xyz_0(0),0,&wrt).run();
  }
  void getRVECoords(RVE * rve, apf::DynamicVector & xyz)
  {
    int ndofs = rve->getDim() * rve->numNodes();
    int sz = xyz.getSize();
    if(sz != ndofs)
      xyz.setSize(ndofs);
    amsi::WriteOp wrt;
    amsi::AccumOp acm;
    amsi::ToArray(rve->getNumbering(),rve->getMesh()->getCoordinateField(),&xyz(0),0,&wrt).run();
    amsi::ToArray(rve->getNumbering(),rve->getUField(),&xyz(0),0,&acm).run();
  }
}
