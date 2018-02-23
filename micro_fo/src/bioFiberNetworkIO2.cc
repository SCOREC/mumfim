#include "bioFiberNetworkIO2.h"
#include "bioFiber.h"
#include "bioFiberNetwork2.h"
#include "bioFiberReactions.h"
#include <apfMeshUtil.h>
#include <cassert>
#include <fstream>
#include <iostream>
namespace bio
{
  class NetworkLoader
  {
  protected:
    apf::Mesh2 * msh;
    apf::MeshTag * prd_tg;
    apf::MeshTag * ord_tg;
    int vrt_cnt;
    int edg_cnt;
    std::vector<apf::MeshEntity*> vrts;
    std::vector<apf::MeshEntity*> edgs;
    void processVertLine(std::istream &);
    void processEdgeLine(std::istream &);
    void processPeriodicity(std::istream &);
    apf::MeshEntity * processVertex(int, const apf::Vector3 &);
  public:
    NetworkLoader()
      : msh(NULL)
      , prd_tg(NULL)
      , ord_tg(NULL)
      , vrt_cnt(0)
      , edg_cnt(0)
      , vrts()
      , edgs()
    {}
    apf::Mesh2 * fromStream(std::istream &);
  };
  apf::Mesh2 * loadFromStream(std::istream & strm)
  {
    return NetworkLoader().fromStream(strm);
  }
  apf::Mesh2 * loadFromFile(const std::string & fnm)
  {
    std::fstream strm(fnm);
    return NetworkLoader().fromStream(strm);
  }
  apf::Mesh2 * NetworkLoader::fromStream(std::istream & is)
  {
    msh = amsi::makeNullMdlEmptyMesh();
    int nn = 0;
    int ne = 0;
    int np = 0;
    is >> nn >> ne >> np;
    if(np)
    {
      prd_tg = msh->createIntTag("periodic_bcs",1);
      assert(sizeof(void*) <= sizeof(int));
    }
    ord_tg = msh->createIntTag("id", 1);
    for(int ii = 0; ii < nn; ++ii)
      processVertLine(is);
    for(int ii = 0; ii < ne; ++ii)
      processEdgeLine(is);
    for(int ii = 0; ii < np; ++ii)
      processPeriodicity(is);
    msh->acceptChanges();
    apf::printStats(msh);
    assert(nn == vrt_cnt);
    assert(ne == edg_cnt);
    return msh;
  }
  void NetworkLoader::processVertLine(std::istream & is)
  {
    apf::Vector3 crd;
    is >> crd[0] >> crd[1] >> crd[2];
    apf::MeshEntity * vrt = msh->createVertex(NULL,crd,apf::Vector3(0,0,0));
    msh->setIntTag(vrt, ord_tg, &vrt_cnt);
    vrt_cnt++;
    vrts.push_back(vrt);
  }
  void NetworkLoader::processEdgeLine(std::istream & is)
  {
    int n0 = -1;
    int n1 = -1;
    is >> n0 >> n1;
    apf::MeshEntity * v[2];
    v[0] = vrts[n0];
    v[1] = vrts[n1];
    apf::MeshEntity * edg = msh->createEntity(apf::Mesh::EDGE, NULL, &v[0]);
    msh->setIntTag(edg, ord_tg, &edg_cnt);
    edg_cnt++;
    edgs.push_back(edg);
  }
  void NetworkLoader::processPeriodicity(std::istream & is)
  {
    int n0 = -1;
    int n1 = -1;
    is >> n0 >> n1;
    msh->setIntTag(vrts[n0-1],prd_tg,reinterpret_cast<int*>(&vrts[n0-1]));
    msh->setIntTag(vrts[n1-1],prd_tg,reinterpret_cast<int*>(&vrts[n1-1]));
  }
}
