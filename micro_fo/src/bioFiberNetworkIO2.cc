#include "bioFiberNetworkIO2.h"
#include "bioFiber.h"
#include "bioFiberNetwork2.h"
#include "bioFiberReactions.h"
#include "bioUtil.h"
#include <cassert>
#include <fstream>
#include <iostream>
namespace bio
{
  FiberNetwork * loadFromFile(const std::string & fnm)
  {
    std::fstream strm(fnm);
    return NetworkLoader().fromStream(strm);
  }
  FiberNetwork * NetworkLoader::fromStream(std::istream & is)
  {
    apf::Mesh2 * msh = makeNullMdlEmptyMesh();
    rct_tg = msh->createIntTag("fiber_reaction",1);
    int et = 0;
    int nr = 0;
    int nn = 0;
    int ne = 0;
    int np = 0;
    is >> et >> nr >> nn >> ne >> np;
    if(np)
    {
      prd_tg = msh->createLongTag("periodic_bcs",1);
      assert(sizeof(void*) <= sizeof(long));
    }
    for(int ii = 0; ii < nr; ++ii)
      processReactionLine(is);
    for(int ii = 0; ii < nn; ++ii)
      processVertLine(is,msh);
    for(int ii = 0; ii < ne; ++ii)
      processEdgeLine(is,msh);
    for(int ii = 0; ii < np; ++ii)
      processPeriodicity(is,msh);
    msh->acceptChanges();
    apf::printStats(msh);
    return new FiberNetwork(msh,(FiberMember)et);
  }
  void NetworkLoader::processReactionLine(std::istream & is)
  {
    int tp = -1;
    is >> tp;
    if(tp == FiberConstitutive::linear)
    {
      LinearReaction * rct = new LinearReaction;
      is >> rct->E;
      rct->fiber_area = 0.0;
      rctns.push_back(rct);
    }
    else if(tp == FiberConstitutive::nonlinear)
    {
      NonlinearReaction * rct = new NonlinearReaction;
      is >> rct->E >> rct->B >> rct->length_ratio_trns;
      rct->fiber_area = 0.0;
      rctns.push_back(rct);
    }
  }
  void NetworkLoader::processVertLine(std::istream & is, apf::Mesh2 * msh)
  {
    apf::Vector3 crd;
    is >> crd[0] >> crd[1] >> crd[2];
    vrts.push_back(msh->createVertex(NULL,crd,apf::Vector3(0,0,0)));
  }
  void NetworkLoader::processEdgeLine(std::istream & is, apf::Mesh2 * msh)
  {
    int n0 = -1;
    int n1 = -1;
    int rct_tp = -1;
    is >> n0 >> n1 >> rct_tp;
    apf::MeshEntity * v[2];
    v[0] = vrts[n0-1];
    v[1] = vrts[n1-1];
    msh->setIntTag(msh->createEntity(apf::Mesh::EDGE,NULL,&v[0]),rct_tg,&rct_tp);
  }
  void NetworkLoader::processPeriodicity(std::istream & is, apf::Mesh2 * msh)
  {
    int n0 = -1;
    int n1 = -1;
    is >> n0 >> n1;
    msh->setLongTag(vrts[n0-1],prd_tg,reinterpret_cast<long*>(vrts[n0-1]));
    msh->setLongTag(vrts[n1-1],prd_tg,reinterpret_cast<long*>(vrts[n1-1]));
  }
}
