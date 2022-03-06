#include "FiberNetwork.h"
#include <iostream>
#include <cassert> // assert
#include <numeric> // accumulate
#include <string>
#include <mth.h>
#include <apfMDS.h>
#include <apfConvert.h>
#include <apfMeshUtil.h> // amsi::makeNullMdlEmptyMesh()
#include <gmi.h>
#include <cstdlib>
#include <memory>
#include "FiberNetworkIO.h"
namespace mumfim
{
  FiberNetworkReactions::FiberNetworkReactions(apf::Mesh2 * msh, std::istream & strm)
  {
    mReactionsList = {};
    loadParamsFromStream(msh, strm, std::back_inserter(mReactionsList));
  }
  FiberNetworkReactions::~FiberNetworkReactions()
  {
    for(auto & reaction : mReactionsList)
    {
      delete reaction;
      reaction = nullptr;
    }
  }
  FiberReaction &  FiberNetworkReactions::operator[](size_t idx)
    {
      return *(mReactionsList[idx]);
    }
  //const FiberReaction & FiberNetworkReactions::operator[](size_t idx) const
  //{
  //  return *(mReactionsList[idx]);
  //}
  FiberNetworkBase::FiberNetworkBase(mesh_ptr_type mesh,
                                     reaction_ptr_type reactions)
      : mMesh(std::move(mesh)), mReactions(std::move(reactions)), mRVEType(0)
  {
    if (mMesh == nullptr || mReactions == nullptr)
    {
      std::cerr << "Attempting to create a fiber network with a null mesh, or "
                   "no reactions!";
      std::cerr << "This is invalid.\n";
      std::exit(EXIT_FAILURE);
    }
    else
    {
    }
  }
  FiberNetworkBase::FiberNetworkBase(const FiberNetworkBase & other) : 
    mMesh(mumfim::make_unique(apf::createMdsMesh(gmi_load(".null") , other.mMesh.get()))),
    mReactions(other.mReactions), mRVEType(other.mRVEType) 
  {
  }
  FiberNetworkBase::~FiberNetworkBase() { };
  FiberNetwork::FiberNetwork(mesh_ptr_type mesh, reaction_ptr_type reactions)
      : FiberNetworkBase(std::move(mesh), std::move(reactions))
      , u(nullptr)
      , xpufnc(nullptr)
      , xpu(nullptr)
      , du(nullptr)
      , udof(nullptr)
      , ucnt(0)
      , tp(FiberMember::truss)
      , scale_factor(1)
  {
    u  = apf::createLagrangeField(mMesh.get(),"u",apf::VECTOR,1);
    du = apf::createLagrangeField(mMesh.get(),"du",apf::VECTOR,1);
    v  = apf::createLagrangeField(mMesh.get(),"v",apf::VECTOR,1);
    a  = apf::createLagrangeField(mMesh.get(),"a",apf::VECTOR,1);
    f  = apf::createLagrangeField(mMesh.get(),"f",apf::VECTOR,1);
    xpufnc = new amsi::XpYFunc(mMesh.get()->getCoordinateField(),u);
    xpu = apf::createUserField(mMesh.get(),"xpu",apf::VECTOR,apf::getShape(u),xpufnc);
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
  FiberNetwork::FiberNetwork(const FiberNetwork& net) : FiberNetworkBase(net)
  {
    assert(net.mMesh);
    u = mMesh->findField(apf::getName(net.u));
    du = mMesh->findField(apf::getName(net.du));
    a = mMesh->findField(apf::getName(net.a));
    v = mMesh->findField(apf::getName(net.v));
    f = mMesh->findField(apf::getName(net.f));
    xpu = mMesh->findField(apf::getName(net.xpu));
    xpufnc = new amsi::XpYFunc(mMesh->getCoordinateField(), u);
    apf::updateUserField(xpu, xpufnc);
    // not a clean way to do this...we can either assume that the naming scheme
    // of the numbering won't change, or that we always have the first numbering
    udof = mMesh->findNumbering(apf::getName(net.udof));
    assert(udof);
    assert(apf::getField(udof));
    tp = net.tp;
    // ordering of reactions should be same in copy
    ucnt = net.ucnt;
    vdof = net.vdof; 
    adof = net.adof;
    fdof = net.fdof;
    scale_factor = net.scale_factor;
  }
  FiberNetwork::~FiberNetwork()
  {
    apf::destroyNumbering(udof);
    udof = nullptr;
    apf::destroyField(xpu);
    xpu = nullptr;
    delete xpufnc;
    xpufnc = nullptr;
    apf::destroyField(u);
    u = nullptr;
    //apf::destroyField(du);
    du = nullptr;
  }
  // TODO this function will need to be updated if we ever use partitioned mesh
  // at the microscale
  void get3DOrientationTensor(mumfim::FiberNetwork* fn, double omega[9])
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
  void get2DOrientationTensor(mumfim::FiberNetwork* fn, double const normal[3],
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
