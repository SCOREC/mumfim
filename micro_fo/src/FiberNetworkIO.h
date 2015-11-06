#ifndef BIO_FIBER_NETWORK_IO_H_
#define BIO_FIBER_NETWORK_IO_H_

#include <iostream>
#include <apf.h>
#include <apfMesh2.h>

namespace bio
{
  class NetworkLoader
  {
  protected:
    std::map<int,apf::MeshEntity*> ldd;

    void parseHeader(std::istream &, int & ne);
    void parseLine(std::istream &, int & n1, int & n2, apf::Vector3 & n1c, apf::Vector3 & n2c);
    void processLine(std::istream &, apf::Mesh2 *);
    apf::MeshEntity * processVertex(apf::Mesh2 *, int, const apf::Vector3 &);
    
  public:
    NetworkLoader() : ldd() {}
    apf::Mesh2 * fromStream(std::istream &);
  };
}

#endif
