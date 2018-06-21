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
    apf::Field * dw;
    apf::Numbering * udof;
    apf::Numbering * wdof;
    int ucnt;
    int wcnt;
    FiberMember tp;
    int dim;
    std::vector<FiberReaction*> rctns;
  public:
    /**
     * Construct a FiberNetwork object.
     * @param f A pointer to a fiber network mesh (contains only vertices and edges)
     *          typically loaded using the NetworkLoader classes
     */
    FiberNetwork(apf::Mesh * f);
    ~FiberNetwork();
    /**
     *  Gives the dimensionality of the managed fiber network
     *  @return the dimensionality of the fiber network (2 or 3)
     */
    int getDim() const { return fn->getDimension(); }
    int getDofCount() const { return ucnt + wcnt; }
    FiberMember getFiberMember()      { return tp;   }
    apf::Mesh * getNetworkMesh()      { return fn;   }
    apf::Field * getUField()          { return u;    }
    apf::Field * getdUField()         { return du;   }
    apf::Field * getXpUField()        { return xpu;  }
    apf::Field * getdWField()         { return dw;   }
    apf::Numbering * getUNumbering()  { return udof; }
    apf::Numbering * getdWNumbering() { return wdof; }
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