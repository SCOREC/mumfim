#include "bioFiberNetwork.h"
#include <iostream>
#include <cassert> // assert
#include <numeric> // accumulate
#include <string>
#include <mth.h>
#include <apfMDS.h>
#include <apfConvert.h>
#include <apfMeshUtil.h> // amsi::makeNullMdlEmptyMesh()
#include <gmi.h>

namespace bio
{
  FiberNetwork::FiberNetwork(apf::Mesh * mesh)
    : fn(mesh)
    , u(NULL)
    , xpufnc(NULL)
    , xpu(NULL)
    , du(NULL)
    , udof(NULL)
    , ucnt(0)
    , tp(FiberMember::truss)
    , rve_tp(0)
    , scale_factor(1)
  {
    assert(mesh);
    du = apf::createLagrangeField(fn,"du",apf::VECTOR,1);
    u  = apf::createLagrangeField(fn,"u",apf::VECTOR,1);
    v  = apf::createLagrangeField(fn,"v",apf::VECTOR,1);
    a  = apf::createLagrangeField(fn,"a",apf::VECTOR,1);
    f  = apf::createLagrangeField(fn,"f",apf::VECTOR,1);
    xpufnc = new amsi::XpYFunc(fn->getCoordinateField(),u);
    xpu = apf::createUserField(fn,"xpu",apf::VECTOR,apf::getShape(u),xpufnc);
    apf::zeroField(du);
    apf::zeroField(u);
    // FIXME...we have to zero the field here or we run into problems copying an
    // empty field for the implicit analysis.
    // we don't zero the velocity and acceleration fields here because
    // we only need them in explicit analysis...only zero when an explicit analysis is created
    apf::zeroField(v);
    apf::zeroField(a);
    apf::zeroField(f);
    // USING 'u' HERE IS VERY IMPORTANT,
    // APPLYDEFORMATIONGRADIENT USES THE
    // FIELD RETRIEVED FROM THIS NUMBERING TO
    // MODIFY, THIS SHOULD BE FIXED BECAUSE
    // IT IS TOO FRAGILE
    udof = apf::createNumbering(u);
    ucnt = apf::NaiveOrder(udof);
    vdof = apf::createNumbering(v);
    adof = apf::createNumbering(a);
    fdof = apf::createNumbering(f);
    apf::NaiveOrder(vdof);
    apf::NaiveOrder(adof);
    apf::NaiveOrder(fdof);
  }
  FiberNetwork::FiberNetwork(const FiberNetwork& net)
  {
    assert(net.fn);
    // note model is deleted with the native mesh
    fn = apf::createMdsMesh(gmi_load(".null") , net.fn);
    u = fn->findField(apf::getName(net.u));
    du = fn->findField(apf::getName(net.du));
    a = fn->findField(apf::getName(net.a));
    v = fn->findField(apf::getName(net.v));
    f = fn->findField(apf::getName(net.f));
    xpu = fn->findField(apf::getName(net.xpu));
    xpufnc = new amsi::XpYFunc(fn->getCoordinateField(), u);
    apf::updateUserField(xpu, xpufnc);
    // not a clean way to do this...we can either assume that the naming scheme
    // of the numbering won't change, or that we always have the first numbering
    udof = fn->findNumbering(apf::getName(net.udof));
    assert(udof);
    assert(apf::getField(udof));
    tp = net.tp;
    rve_tp = net.rve_tp;
    // ordering of reactions should be same in copy
    rctns.resize(net.rctns.size());
    std::copy(net.rctns.begin(), net.rctns.end(), rctns.begin());
    ucnt = net.ucnt;
    vdof = net.vdof; 
    adof = net.adof;
    fdof = net.fdof;
    scale_factor = net.scale_factor;
  }
  FiberNetwork::~FiberNetwork()
  {
    apf::destroyNumbering(udof);
    udof = NULL;
    apf::destroyField(xpu);
    xpu = NULL;
    delete xpufnc;
    xpufnc = NULL;
    apf::destroyField(u);
    u = NULL;
    //apf::destroyField(du);
    du = NULL;
    fn->destroyNative();
    apf::destroyMesh(fn);
    fn = NULL;
  }
  // TODO this function will need to be updated if we ever use partitioned mesh
  // at the microscale
  void get3DOrientationTensor(bio::FiberNetwork* fn, double omega[9])
  {
    for (int i = 0; i < 9; ++i) {
      omega[i] = 0;
    }
    apf::Mesh* network = fn->getNetworkMesh();
    unsigned int edgeCount = 0;
    // loop over all of the edges in the network
    apf::MeshIterator* it = network->begin(1);
    while (apf::MeshEntity* edge = network->iterate(it)) {
      // edges only have two adjacent vertices
      apf::MeshEntity* adjVerts[2];
      network->getDownward(edge, 0, adjVerts);
      // get vector between two verts of edge
      apf::Vector3 p0, p1;
      apf::getVector(fn->getXpUField(), adjVerts[0], 0, p0);
      apf::getVector(fn->getXpUField(), adjVerts[1], 0, p1);
      apf::Vector3 edgeUnit = (p1 - p0).normalize();
      apf::Matrix3x3 temp = apf::tensorProduct(edgeUnit, edgeUnit);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          omega[i * 3 + j] = omega[i * 3 + j] + temp[i][j];
        }
      }
      ++edgeCount;
    }
    network->end(it);
    // normalize r by the number of edges
    for (int i = 0; i < 9; ++i) {
      omega[i] = omega[i] / edgeCount;
    }
  }
  // TODO this function will need to be updated if we ever use partitioned mesh
  // at the microscale
  void get2DOrientationTensor(bio::FiberNetwork* fn, double const normal[3],
                              double omega[9])
  {
    for (int i = 0; i < 9; ++i) {
      omega[i] = 0;
    }
    apf::Mesh* network = fn->getNetworkMesh();
    // not super efficent, but then we don't have to redefine dot product and
    // multiplications
    apf::Vector3 n(normal);
    // normalize the vector defining the plane because it not gaurenteed to be
    // normalized from simmodeler
    // FIXME it may be better to do this at macroscale because then it is only
    // normalized once...
    n = n.normalize();
    unsigned int edgeCount = 0;
    apf::MeshIterator* it = network->begin(1);
    while (apf::MeshEntity* edge = network->iterate(it)) {
      apf::MeshEntity* adjVerts[2];
      network->getDownward(edge, 0, adjVerts);
      apf::Vector3 p0, p1;
      apf::getVector(fn->getXpUField(), adjVerts[0], 0, p0);
      apf::getVector(fn->getXpUField(), adjVerts[1], 0, p1);
      apf::Vector3 edgeUnit = (p1 - p0).normalize();
      // project edge unit into plane defined by normal
      // e.g. t = edgeUnit-(edgeUnit.n)n (subtract out the normal part)
      // we could use apf::reject here in future
      apf::Vector3 t = edgeUnit - n * (edgeUnit * n);
      t = t.normalize();
      apf::Matrix3x3 temp = apf::tensorProduct(t, t);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          omega[i * 3 + j] = omega[i * 3 + j] + temp[i][j];
        }
      }
      ++edgeCount;
    }
    network->end(it);
    for (int i = 0; i < 9; ++i) {
      omega[i] = omega[i] / edgeCount;
    }
  }
}
