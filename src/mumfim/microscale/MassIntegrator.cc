#include "MassIntegrator.h"
#include <cassert>
namespace mumfim
{
  void MassIntegrator::inElement(apf::MeshElement * me)
  {
    mass_elmt = apf::createElement(mass_matrix_fld, me);
    nnd = apf::countNodes(mass_elmt);
    cmp_per_nd = apf::countComponents(mass_matrix_fld);
    nedofs = nnd * cmp_per_nd;
    mass_elem = apf::DynamicMatrix(nedofs, nedofs);
    mass_elem.zero();
    mesh->getDownward(apf::getMeshEntity(me), 0, &nodes[0]);
  }
  // integrate density to get mass
  void MassIntegrator::atPoint(apf::Vector3 const & p, double w, double dV)
  {
    // we don't currently support not lumping the mass matrix
    // becuase we use the mass matrix as a field in the element by element
    // explicit method
    assert(massLumpType != MassLumpType::None);
    apf::NewArray<double> N;
    apf::getShapeValues(mass_elmt, p, N);
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
          mass_elem(dof1, dof2) += density * N1 * N2 * w * fiber_area * dV;
        }
      }
    }
  }
  // assemble mass matrix
  void MassIntegrator::outElement()
  {
    for (unsigned int i = 0; i < nnd; ++i)
    {
      double val = apf::getScalar(mass_matrix_fld, nodes[i], 0);
      int dof1 = i * cmp_per_nd;
      apf::setScalar(mass_matrix_fld, nodes[i], 0, val + mass_elem(dof1, dof1));
    }
    apf::destroyElement(mass_elmt);
  }
}  // namespace mumfim
