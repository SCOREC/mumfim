#include "bioFiberNetwork.h"
#include <iostream>
#include <cassert> // assert
#include <numeric> // accumulate
#include <string>
#include <mth.h>
namespace bio
{
  FiberNetwork::FiberNetwork(apf::Mesh * f)
    : fn(f)
    , u(NULL)
    , xpufnc(NULL)
    , xpu(NULL)
    , du(NULL)
    , udof(NULL)
    , ucnt(0)
    , tp(FiberMember::truss)
  {
    assert(f);
    du = apf::createLagrangeField(fn,"du",apf::VECTOR,1);
    u  = apf::createLagrangeField(fn,"u",apf::VECTOR,1);
    xpufnc = new amsi::XpYFunc(fn->getCoordinateField(),u);
    xpu = apf::createUserField(fn,"xpu",apf::VECTOR,apf::getShape(u),xpufnc);
    apf::zeroField(du);
    apf::zeroField(u);
    // USING 'u' HERE IS VERY IMPORTANT,
    // APPLYDEFORMATIONGRADIENT USES THE
    // FIELD RETRIEVED FROM THIS NUMBERING TO
    // MODIFY, THIS SHOULD BE FIXED BECAUSE
    // IT IS TOO FRAGILE
    udof = apf::createNumbering(u);
    ucnt = apf::NaiveOrder(udof);
  }
  FiberNetwork::~FiberNetwork()
  {
    apf::destroyNumbering(udof);
    apf::destroyField(xpu);
    delete xpufnc;
    apf::destroyField(u);
    apf::destroyField(du);
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
