#ifndef BIO_LINEARTISSUE_H_
#define BIO_LINEARTISSUE_H_
#include <apfSimFEA.h>
namespace bio
{
  class LinearTissue : public amsi::apfSimFEA
  {
  public:
    LinearTissue(pGModel imdl, pParMesh imsh, pACase ipd, pACase iss, MPI_Comm cm);
    virtual void UpdateDOFs(const double * sol);
    apf::Field * getField();
  };
}
#endif
