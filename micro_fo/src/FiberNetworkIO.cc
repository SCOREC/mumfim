#include "FiberNetworkIO.h"
#include "apfUtil.h"
#include "FiberNetwork.h"
#include <cassert>
#include <fstream>
#include <iostream>
namespace bio
{
  FiberNetwork * loadFromFile(const std::string & fnm)
  {
    std::fstream strm(fnm);
    return new FiberNetwork(NetworkLoader().fromStream(strm));
  }
  
  apf::Mesh2 * NetworkLoader::fromStream(std::istream & is)
  {
    apf::Mesh2 * msh = makeNullMdlEmptyMesh();
    int num_edges = 0;
    parseHeader(is,num_edges);
    for(int ii = 0; ii < num_edges; ii++)
      processLine(is,msh);
    msh->acceptChanges();
    apf::printStats(msh);
    return msh;
  }

  void NetworkLoader::parseHeader(std::istream & is, int & num_edges)
  {
    int nn, ndof;
    is >> nn >> ndof >> num_edges;
  }

  void NetworkLoader::parseLine(std::istream & is,
				int & n1, int & n2,
				apf::Vector3 & n1p, apf::Vector3 & n2p)
  {
    int edge_id = -1;
    double n1d[3];
    double n2d[3];
    is >> edge_id >> n1 >> n2 >> n1d[0] >> n1d[1] >> n1d[2] >> n2d[0] >> n2d[1] >> n2d[2];
    n1p.fromArray(n1d);
    n2p.fromArray(n2d);
  }

  void NetworkLoader::processLine(std::istream & is, apf::Mesh2 * msh)
  {
    int n1 = -1;
    int n2 = -1;
    apf::Vector3 n1c;
    apf::Vector3 n2c;
    parseLine(is,n1,n2,n1c,n2c);
    apf::MeshEntity * v[2];
    v[0] = processVertex(msh,n1,n1c);
    v[1] = processVertex(msh,n2,n2c);
    msh->createEntity(apf::Mesh::EDGE,NULL,&v[0]);
  }

  apf::MeshEntity * NetworkLoader::processVertex(apf::Mesh2 * msh,
						 int id,
						 const apf::Vector3 & coord)
  {
    apf::MeshEntity * result = NULL;
    if(ldd.count(id))
      result = ldd[id];
    else
    {
      result = msh->createVertex(NULL,coord,apf::Vector3(0,0,0));
      ldd[id] = result;
    }
    return result;
  }





  

}
