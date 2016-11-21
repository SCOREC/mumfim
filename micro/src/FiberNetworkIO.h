#ifndef BIO_FIBER_NETWORK_IO_H_
#define BIO_FIBER_NETWORK_IO_H_

#include <iostream>
#include <apf.h>
#include <apfMesh2.h>

namespace bio
{

  /**
   * A utility class typically used as an anonymous instance simply to process the
   *  loading of the fiber network mesh from different sources.
   * For example in load.cc : \snippet test/load.cc load from file stream
   * @note When additional sources are being added, reevaluate whether this class should be
   *       made into an abstract base class with different derived classes for different sources
   *       since there is only a single source of fiber networks at the moment the decision
   *       is moot at this point in time.
   */
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
