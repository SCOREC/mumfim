#include "bioFiberNetwork.h"
#include <apfDynamicMatrix.h>
namespace bio
{
  template <typename O>
    void getBoundaryVerts(const RVE * rve, apf::Mesh * fn_msh, RVE::side sd, O nds)
  {
    apf::MeshEntity * me = NULL;
    apf::MeshIterator * it = NULL;
    for(it = fn_msh->begin(0); (me = fn_msh->iterate(it)) ; )
    {
      apf::Vector3 crd;
      fn_msh->getPoint(me,0,crd);
      if(rve->onBoundary(crd,sd))
        *nds++ = me;
    }
    fn_msh->end(it);
  }
}
