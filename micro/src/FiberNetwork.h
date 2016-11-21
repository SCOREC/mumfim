#ifndef BIO_FIBER_NETWORK_H_
#define BIO_FIBER_NETWORK_H_
#include "apfUtil.h"
#include "lasCSR.h"
#include "ElementalSystem.h"
#include "FiberReactions.h"
#include "RVE.h"
#include "lasSparskit.h"
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfField.h>
#include <cassert>
namespace bio
{
  /**
   * Responsible for managing the internal state of a single fiber-network
   *  quasistatics simulation. 
   */
  class FiberNetwork
  {
  private:
    apf::Mesh * fn;
    apf::Field * fn_u;
    apf::Field * fn_du;
    apf::Numbering * fn_dof;
    int dim;
  protected:
  public:

    /**
     * Construct a FiberRVE object.
     * @param f A pointer to a fiber network mesh (contains only vertices and edges)
     *          typically loaded using the NetworkLoader classes 
     */
    FiberNetwork(apf::Mesh * f);

    /**
     *  Gives the dimensionality of the managed fiber network
     *  @return the dimensionality of the fiber network (2 or 3)
     */
    int getDim() const { return fn->getDimension(); }

    /**
     * Calculate the current length of each fiber in the network
     * @param lngths A vector containing the lengths of all fibers in the network;
     * @note Possibly implement local caching of the answer instead of recalculating
     */
    void calcFiberLengths(std::vector<double> & lngths) const
    {
      calcDimMeasures(fn,1,lngths);
    }
    apf::Mesh * getNetworkMesh() { return fn; }
    apf::Field * getDisplacementField() { return fn_u; }
    apf::Field * getIncrementalDispField() { return fn_du; }
    apf::Numbering * getNumbering() { return fn_dof; }
  };

  /**
   * TODO (m) Bill : move to FiberNetworkIO files
   */
  FiberNetwork * loadFromFile(const std::string & fnm);

  void assembleElementalSystem(las::skMat * k,
			                         las::skVec * f,
			                         const ElementalSystem * es,
			                         apf::NewArray<int> & dofs);
}
#endif
