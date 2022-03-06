#ifndef MUMFIM_MASS_INTEGRATOR_H_
#define MUMFIM_MASS_INTEGRATOR_H_
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfNumbering.h>
//#include "bioFiberReactions.h"
namespace mumfim
{
  /*
   * Currently only the row sum technique is implemented (referenced in
   * Belytschko, but there is discussion of a lobatto/optimal lumping type in
   * the cook textbook.
   */
  enum class MassLumpType
  {
    None,
    RowSum
  };
  /*
   * compute the consistent mass matrix
   */
  class MassIntegrator : public apf::Integrator
  {
    protected:
    apf::Element * mass_elmt;
    apf::DynamicMatrix mass_elem;
    unsigned int nnd;
    unsigned int cmp_per_nd;
    unsigned int nedofs;
    MassLumpType massLumpType;
    double density;
    double fiber_area;
    apf::Mesh * mesh;
    apf::Field * mass_matrix_fld;
    apf::MeshEntity * nodes[2];
    public:
    MassIntegrator(apf::Field * mass_matrix_fld, double density,
                   double fiber_area, int order,
                   MassLumpType mlt = MassLumpType::None)
        : apf::Integrator(order)
        , massLumpType(mlt)
        , density(density)
        , fiber_area(fiber_area)
        , mesh(apf::getMesh(mass_matrix_fld))
        , mass_matrix_fld(mass_matrix_fld)
    {
    }
    void inElement(apf::MeshElement * me) override;
    void atPoint(apf::Vector3 const & p, double w, double dV) override;
    void outElement() override;
  };
}  // namespace mumfim
#endif
