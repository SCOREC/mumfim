#include "bioFiberNetwork.h"
#include <apfDynamicMatrix.h>
namespace bio
{
  template <typename I, typename O>
    void getBoundaryVerts(const RVE * rve,
                          apf::Mesh * msh,
                          I bgn_vrts,
                          I end_vrts,
                          RVE::side sd,
                          O nds)
  {
    for(auto vrt = bgn_vrts; vrt != end_vrts; ++vrt)
    {
      apf::Vector3 crd;
      msh->getPoint(*vrt,0,crd);
      if(rve->onBoundary(crd,sd))
        *nds++ = *vrt;
    }
  }
}
