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
  FiberNetwork * loadFromFileAndParams(const std::string & fnm)
  {
    std::fstream strm(fnm);
    NetworkLoader ldr;
    FiberNetwork * fn = ldr.fromStream(strm);
    strm.close();
    std::string prms_fnm = fnm + std::string(".params");
    std::fstream prm_strm(prms_fnm);
    ldr.paramsFromStream(prm_strm,std::back_inserter(fn->getFiberReactions()));
    return fn;
  }
  FiberNetwork * NetworkLoader::fromStream(std::istream & is)
  {
    msh = makeNullMdlEmptyMesh();
    int nn = 0;
    int ne = 0;
    int np = 0;
    is >> nn >> ne >> np;
    if(np)
    {
      prd_tg = msh->createLongTag("periodic_bcs",1);
      assert(sizeof(void*) <= sizeof(long));
    }
    for(int ii = 0; ii < nn; ++ii)
      processVertLine(is);
    for(int ii = 0; ii < ne; ++ii)
      processEdgeLine(is);
    for(int ii = 0; ii < np; ++ii)
      processPeriodicity(is);
    msh->acceptChanges();
    apf::printStats(msh);
    return new FiberNetwork(msh);
  }
  void NetworkLoader::processVertLine(std::istream & is)
  {
    apf::Vector3 crd;
    is >> crd[0] >> crd[1] >> crd[2];
    vrts.push_back(msh->createVertex(NULL,crd,apf::Vector3(0,0,0)));
  }
  void NetworkLoader::processEdgeLine(std::istream & is)
  {
    int n0 = -1;
    int n1 = -1;
    is >> n0 >> n1;
    apf::MeshEntity * v[2];
    v[0] = vrts[n0];
    v[1] = vrts[n1];
    edgs.push_back(msh->createEntity(apf::Mesh::EDGE,NULL,&v[0]));
  }
  void NetworkLoader::processPeriodicity(std::istream & is)
  {
    int n0 = -1;
    int n1 = -1;
    is >> n0 >> n1;
    msh->setLongTag(vrts[n0-1],prd_tg,reinterpret_cast<long*>(vrts[n0-1]));
    msh->setLongTag(vrts[n1-1],prd_tg,reinterpret_cast<long*>(vrts[n1-1]));
  }
  void NetworkLoader::processEdgeReactionLine(std::istream & is, int ii)
  {
    int r = -1;
    is >> r;
    msh->setIntTag(edgs[ii],rct_tg,&r);
  }
}
