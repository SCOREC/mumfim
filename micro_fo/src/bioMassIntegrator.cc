#include "bioMassIntegrator.h"
#include <cassert>
namespace bio
{
  void MassIntegrator::inElement(apf::MeshElement * me)
  {
    prim_elmt = apf::createElement(primary_fld, me);
    nnd = apf::countNodes(prim_elmt);
    cmp_per_nd = apf::countComponents(primary_fld);
    nedofs = nnd * cmp_per_nd;
    mass_elem = apf::DynamicMatrix(nedofs, nedofs);
    mass_elem.zero();
    apf::getElementNumbers(nm, apf::getMeshEntity(me), dofs);
    int tg = -1;
    mesh->getIntTag(apf::getMeshEntity(me), rct_tg, &tg);
    fr = frs[tg];
  }
  // integrate density to get mass
  void MassIntegrator::atPoint(apf::Vector3 const & p, double w, double dV)
  {
    apf::NewArray<double> N;
    apf::getShapeValues(prim_elmt, p, N);
    // Want to make sure there isn't an entry
    // for each dof at each node, because that
    // will make the indexing wrong
    assert(N.size() == nnd);
    for (unsigned int j = 0; j < cmp_per_nd; ++j)
    {
      for (unsigned int i = 0; i < nnd; ++i)
      {
        double N1 = N[i];
        int dof1 = i * cmp_per_nd + j;
        for (unsigned int k = 0; k < nnd; ++k)
        {
          double N2 = N[k];
          int dof2 = k * cmp_per_nd + j;
          if (massLumpType == MassLumpType::RowSum) dof2 = dof1;
          // mass = integral(density*Na*Nb*dV)
          mass_elem(dof1, dof2) += fr->fiber_density * N1 * N2 * w * dV;
        }
      }
    }
  }
  // assemble mass matrix
  void MassIntegrator::outElement()
  {
    // set mass for this element number
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    // mass matrix
    // TODO if mass matrix lumped, only assemble into the diagonals
    ops->assemble(mass, nedofs, &dofs[0], nedofs, &dofs[0], &mass_elem(0, 0));
    apf::destroyElement(prim_elmt);
  }
}  // namespace bio
