#ifndef BIO_FIBER_NETWORK_IO2_H_
#define BIO_FIBER_NETWORK_IO2_H_
#include <bioFiberReactions.h>
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
    apf::MeshTag * rct_tg;
    apf::MeshTag * prd_tg;
    apf::Field * mmts;
    std::vector<apf::MeshEntity*> vrts;
    std::vector<FiberReaction*> rctns;
    void processReactionLine(std::istream &);
    void processVertLine(std::istream &, int, apf::Mesh2 *);
    void processEdgeLine(std::istream &, apf::Mesh2 *);
    void processPeriodicity(std::istream &, apf::Mesh2 *);
    apf::MeshEntity * processVertex(apf::Mesh2 *, int, const apf::Vector3 &);
  public:
    NetworkLoader()
      : rct_tg()
      , prd_tg()
      , mmts()
      , vrts()
      , rctns()
    {}
    apf::Mesh2 * fromStream(std::istream &);
  };
}
#endif
