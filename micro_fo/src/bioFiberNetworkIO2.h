#ifndef BIO_FIBER_NETWORK_IO2_H_
#define BIO_FIBER_NETWORK_IO2_H_
#include "bioFiberNetwork.h"
#include "bioFiberReactions.h"
#include <iostream>
#include <apf.h>
#include <apfMesh2.h>
namespace bio
{
  FiberNetwork * loadFromFile(const std::string & fnm);
  FiberNetwork * loadFromFileAndParams(const std::string & fnm);
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
    apf::Mesh2 * msh;
    apf::MeshTag * rct_tg;
    apf::MeshTag * prd_tg;
    std::vector<apf::MeshEntity*> vrts;
    std::vector<apf::MeshEntity*> edgs;
    void processVertLine(std::istream &);
    void processEdgeLine(std::istream &);
    void processPeriodicity(std::istream &);
    apf::MeshEntity * processVertex(int, const apf::Vector3 &);
    template <typename O>
      void processReactionLine(std::istream &, O);
    void processEdgeReactionLine(std::istream &,int);
  public:
    NetworkLoader()
      : msh(NULL)
      , rct_tg(NULL)
      , prd_tg(NULL)
      , vrts()
      , edgs()
    {}
    FiberNetwork * fromStream(std::istream &);
    template <typename O>
      void paramsFromStream(std::istream &, O);
  };
}
#include "bioFiberNetworkIO2_impl.h"
#endif
