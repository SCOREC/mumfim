#include "bioFiberNetwork.h"
#include <apfDynamicMatrix.h>
namespace bio
{
  template <typename O>
    void getBoundaryVerts(RVE * rve, FiberNetwork * fn, RVE::side sd, O nds)
  {
    apf::Mesh * fn_msh = fn->getNetworkMesh();
    apf::Field * u = fn->getUField();
    apf::FieldShape * shp = apf::getShape(u);
    apf::MeshEntity * me = NULL;
    apf::MeshIterator * it = NULL;
    for(it = fn_msh->begin(0); me = fn_msh->iterate(it);)
    {
      int cnt = shp->countNodesOn(fn_msh->getType(me));
      for(int nd = 0; nd < cnt; nd++)
      {
        apf::Vector3 crd;
        fn_msh->getPoint(me,nd,crd);
        if(rve->onBoundary(crd,sd))
          *nds++ = me;
      }
    }
    fn_msh->end(it);
  }
}
