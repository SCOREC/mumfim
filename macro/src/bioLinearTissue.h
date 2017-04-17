#ifndef BIO_LINEARTISSUE_H_
#define BIO_LINEARTISSUE_H_
#include <apfSimFEA.h>
namespace bio
{
  class LinearTissue : public amsi::apfSimFEA
  {
  private:
    apf::Field * stress_ip_field;
    apf::Field * strain_ip_field;
    apf::Integrator * strss_strn_itgr;
  public:
    LinearTissue(pGModel imdl, pParMesh imsh, pACase ipd, apf::Field * strs, apf::Field * strn, MPI_Comm cm);
    virtual void UpdateDOFs(const double * sol);
    virtual void Assemble(amsi::LAS * las);
    apf::Field * getField();
  };
}
#endif
