#include "FiberNetwork.h"
#include "FiberReactions.h"
#include "las.h"
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
    fn_u = apf::createLagrangeField(f,"u",apf::VECTOR,1);
    fn_du = apf::createLagrangeField(f,"du",apf::VECTOR,1);
    fn_dof = apf::createNumbering(fn_du);
  }
  void assembleElementalSystem(las::skMat * k,
			                         las::skVec * f,
			                         const ElementalSystem * es,
			                         apf::NewArray<int> & dofs)
  {
    int nedofs = es->nedofs();
    las::setVecValues(f,es->getfe(),dofs,nedofs,true);
    las::setMatValues(k,es->getKe(),dofs,nedofs,dofs,nedofs,true);
  }
}
