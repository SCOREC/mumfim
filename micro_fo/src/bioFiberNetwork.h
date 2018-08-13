#ifndef BIO_FIBER_NETWORK_H_
#define BIO_FIBER_NETWORK_H_
#include "bioFiberReactions.h"
#include "bioFiber.h"
#include <apfFunctions.h> // amsi
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <string>
#include <vector>
namespace bio
{
  struct FiberNetworkReactions
  {
    apf::Mesh2 * msh;
    std::vector<FiberReaction*> rctns;
    std::string fileName;
  };
  /**
   * Responsible for managing the internal state of a single fiber-network
   *  quasistatics simulation.
   */
  class FiberNetwork
  {
  protected:
    apf::Mesh * fn;
    apf::Field * u;
    amsi::XpYFunc * xpufnc;
    apf::Field * xpu;
    apf::Field * du;
    apf::Numbering * udof;
    int ucnt;
    FiberMember tp;
    std::vector<FiberReaction*> rctns;
    // the rve filename map
    int rve_tp;
  public:
    /**
     * Construct a FiberNetwork object.
     * @param f A pointer to a fiber network mesh (contains only vertices and edges)
     *          typically loaded using the NetworkLoader classes
     */
    FiberNetwork(apf::Mesh * f);
    FiberNetwork(const FiberNetwork & fn);
    ~FiberNetwork();
    /**
     *  Gives the dimensionality of the managed fiber network
     *  @return the dimensionality of the fiber network (2 or 3)
     */
    int getDim() const { return fn->getDimension(); }
    int getDofCount() const { return ucnt; }
    /*
     * returns the type of rve (e.g. filename as an integer)
     * the mapping between this integer and the filename can
     * be found in the rve_tp log
     */
    int getRVEType()                  { return rve_tp;}
    void setRVEType(int rve_t)        { rve_tp = rve_t;}
    FiberMember getFiberMember()      { return tp;   }
    apf::Mesh * getNetworkMesh()      { return fn;   }
    apf::Field * getUField()          { return u;    }
    apf::Field * getdUField()         { return du;   }
    apf::Field * getXpUField()        { return xpu;  }
    apf::Numbering * getUNumbering()  { return udof; }
    std::vector<FiberReaction*> & getFiberReactions() { return rctns; }
  };

  /**
   * get the mean orientation vector of the rve
   * \param network pointer to the the fiber network
   * \param \out omega the orientation tensor
   */
  void get3DOrientationTensor(bio::FiberNetwork * network, double omega[9]);
  /**
   * get the mean orientation vector of the rve projected into the plane defined by the normal
   * \param network pointer to the the fiber network
   * \param normal the normal vector to the plane to project the vectors into
   * \param \out omega orientation tensor
   */
  void get2DOrientationTensor(bio::FiberNetwork* network, double const normal[3], double omega[9]);
}
#endif
