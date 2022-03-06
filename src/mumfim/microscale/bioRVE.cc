#include "bioRVE.h"
#include "bioFiberNetwork.h"
#include "lasSparskit.h"
#include <apfFunctions.h>
#include <apfMeshIterator.h>
#include <apfMeshUtil.h>
#include <apfMDS.h>
#include <apfConvert.h>
namespace mumfim
{
  template <typename O>
  static void originCenterCube(O cbe_crnrs, double crd)
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
  static void originCenterSquare(O cbe_crnrs, double crd)
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
  , cbe_dof_cnt(0)
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
    cbe_dof_cnt = apf::NaiveOrder(cbe_dof);
    assert(dim == 3 ? cbe_dof_cnt == 24 : cbe_dof_cnt == 8);
    cbe_u_e = apf::createElement(cbe_u,cbe_e);
  }
  RVE::RVE(const RVE & rve) {
    dim = rve.dim;
    crd = rve.crd;
    cbe = amsi::makeNullMdlEmptyMesh();
    apf::convert(rve.cbe, static_cast<apf::Mesh2*>(cbe));
    apf::MeshIterator * it = cbe->begin(3);
    cbe_e = cbe->iterate(it);
    cbe->end(it);
    cbe_u = cbe->findField(apf::getName(rve.cbe_u));
    cbe_u_e = apf::createElement(cbe_u,cbe_e);
    cbe_du = cbe->findField(apf::getName(rve.cbe_du));
    cbe_xpu = cbe->findField(apf::getName(rve.cbe_xpu));
    xpufnc = new amsi::XpYFunc(cbe->getCoordinateField(), cbe_u);
    apf::updateUserField(cbe_xpu, xpufnc);
    cbe_dof = cbe->findNumbering(apf::getName(rve.cbe_dof));
    assert(cbe_dof);
    assert(apf::getField(cbe_dof));
  }
  RVE::~RVE()
  {
    apf::destroyElement(cbe_u_e);
    cbe_u_e = NULL;
    apf::destroyField(cbe_u);
    cbe_u = NULL;
    apf::destroyField(cbe_du);
    cbe_du = NULL;
    apf::destroyField(cbe_xpu);
    cbe_xpu = NULL;
    apf::destroyNumbering(cbe_dof);
    delete xpufnc;
    xpufnc = NULL;
    cbe->destroyNative();
    apf::destroyMesh(cbe);
    cbe = NULL;
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
}
