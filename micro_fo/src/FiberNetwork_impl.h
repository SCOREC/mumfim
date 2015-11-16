#include "FiberNetwork_impl.h"

namespace bio
{
  template <typename T, typename OP>
    void assembleLinearSystem(const apf::Mesh * msh,
			      apf::Numbering * num,
			      int dim,
			      T * es_gn,
			      LinearSystem * ls)
  {
    int nedofs = 0;
    apf::MeshEntity * edge = NULL;
    apf::MeshIterator * it = NULL;
    for(it = msh->begin(dim); me = msh->iterate(it);)
    {
      ElementalSystem * es = es_gn->OP(apf::createMeshElement(msh,me));
      apf::NewArray<int> dofs;
      apf::getElementNumbers(num,me,dofs);
      assembleElementalSystem(es,dofs,ls);
    }
    msh->end(it);
  }
}
