#include "bioRVE2.h"
#include "bioFiberNetwork.h"
#include "lasSparskit.h"
#include <apfFunctions.h>
#include <apfMeshUtil.h>
namespace bio
{
  template <typename O>
  void originCenterCube(O cbe_crnrs, double crd)
  {
    *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd * -1.0,crd * -1.0);
    *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd * -1.0,crd * -1.0);
    *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd * -1.0,crd *  1.0);
    *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd * -1.0,crd *  1.0);
    *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd *  1.0,crd * -1.0);
    *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd *  1.0,crd * -1.0);
    *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd *  1.0,crd *  1.0);
    *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd *  1.0,crd *  1.0);
  }
  template <typename O>
  void originCenterSquare(O cbe_crnrs, double crd)
  {
    *cbe_crnrs++ = apf::Vector3(-crd,-crd, 0.0);
    *cbe_crnrs++ = apf::Vector3( crd,-crd, 0.0);
    *cbe_crnrs++ = apf::Vector3( crd, crd, 0.0);
    *cbe_crnrs++ = apf::Vector3(-crd, crd, 0.0);
  }
  RVE::RVE(int d)
  : dim(d)
  , crd(0.5)
  , cbe(NULL)
  , cbe_e(NULL)
  , cbe_u_e(NULL)
  , cbe_u(NULL)
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
    cbe_dof = apf::createNumbering(cbe_u);
    apf::NaiveOrder(cbe_dof);
    cbe_u_e = apf::createElement(cbe_u,cbe_e);
  }
  void forwardRVEDisplacement(RVE * rve, FiberNetwork * fn)
  {
    apf::Element * cbe_e = rve->getElement();
    apf::Field * fn_u = fn->getUField();
    InterpOnto forward_u(cbe_e,fn_u);
    forward_u.apply(fn_u);
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
    ApplySolution(rve->getNumbering(),&du[0],0,true).apply(rve->getField());
  }
  void getRVEDisplacement(RVE * rve, apf::DynamicVector & u)
  {
    int ndofs = rve->getDim() * rve->numNodes();
    int sz = u.getSize();
    if(sz != ndofs)
      u.setSize(ndofs);
    amsi::WriteOp wrt;
    amsi::ToArray(rve->getNumbering(),rve->getField(),&u(0),0,&wrt).run();
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
    amsi::ToArray(rve->getNumbering(),rve->getField(),&xyz(0),0,&acm).run();
  }
}
