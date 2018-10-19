#ifndef BIO_MASS_INTEGRATOR_H_
#define BIO_MASS_INTEGRATOR_H_
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfNumbering.h>
#include <las.h>
#include <lasConfig.h>
namespace bio
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
    // assume this is a step field (one value per element)
    apf::Field * density_fld;
    // we will use the same field interpolations as the primary field
    // for a mechanics analysis this is the displacements/velocities
    apf::Field * primary_fld;
    apf::Numbering * nm;
    apf::Element * dens_elmt;
    apf::Element * prim_elmt;
    apf::DynamicMatrix mass_elem;
    unsigned int nnd;
    unsigned int cmp_per_nd;
    unsigned int nedofs;
    apf::NewArray<int> dofs;
    las::Mat * mass;
    MassLumpType massLumpType;
    public:
    MassIntegrator(apf::Numbering * nm,
                   apf::Field * density,
                   apf::Field * primary,
                   las::Mat * mass,
                   int order,
                   MassLumpType mlt = MassLumpType::None)
        : apf::Integrator(order)
        , density_fld(density)
        , primary_fld(primary)
        , nm(nm)
        , mass(mass)
        , massLumpType(mlt)
    {
    }
    void inElement(apf::MeshElement * me) override;
    void atPoint(apf::Vector3 const & p, double w, double dV) override;
    void outElement() override;
  };
}  // namespace bio
#endif
