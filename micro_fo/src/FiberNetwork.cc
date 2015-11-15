#include "FiberNetwork.h"
#include "TrussIntegrator.h"

#include <apfShape.h>
#include <apfMesh.h>
#include <apfNumbering.h>

#include <cassert>
#include <cmath>
#include <numeric>

namespace bio
{
  FiberNetwork::FiberNetwork(apf::Mesh * f)
    : fn(f)
    , fn_u(NULL)
    , fn_dof(NULL)
    , dim(f->getDimension())
  {
    assert(f);
    fn_u = apf::createLagrangeField(f,"displacement",apf::VECTOR,1);
    fn_dof = apf::createNumbering(fn_u);
  }

  void FiberNetwork::assembleJacobian()
  {
    TrussIntegrator elmnt_stm;

    apf::MeshEntity * edge = NULL;
    apf::MeshIterator * it = NULL;
    for(it = mesh->begin(1); me = mesh->iterate(it);)
    {
      elmnt_stm.process(me);
      apf::DynamicMatrix & ke = elmnt_stm.getKe();
      apf::DynamicMatrix & fe = elmnt_stm.getFe();
    }
    fn->end(it);
  }
}
