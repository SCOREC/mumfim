#include "FiberNetwork.h"
#include "FiberReactions.h"
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

  void FiberNetwork::applySolution(double * sol, int ndof)
  {

  }

  const ElementalSystem * FiberNetwork::calcElementalSystem(apf::MeshElement * melmt)
  {
    // retrieve the integrator responsible for generating the elemental system
    //  for this fiber, use it and return the resulting elemental system
    static TrussIntegrator intgrtr(1,fn_u,new LinearReaction);
    intgrtr.process(melmt);
    return intgrtr.getElementalSystem();
  }

  void assembleElementalSystem(const ElementalSystem * es,
			       const apf::NewArray<int> & dofs,
			       LinearSystem * ls)
  {
    int nedofs = es->nedofs();
    LSOps * ops = ls->getOps();
    ops->addToVector(ls->getVector(),es->getfe(),dofs,nedofs);
    ops->addToMatrix(ls->getMatrix(),es->getKe(),dofs,nedofs);
  }

  /*
  assembleLinearSystem<FiberNetwork,FiberNetwork::calcElementalSystem>(fn->getMesh(),fn,ls);
  assembleLinearSystem<ElementalSystemIntegrator,
		       ElementalSystemIntegrator::calcElementalSystem>(fn->getMesh(),new TrussIntegrator,ls);
  */
}
