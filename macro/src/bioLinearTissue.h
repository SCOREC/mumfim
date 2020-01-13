#ifndef BIO_LINEARTISSUE_H_
#define BIO_LINEARTISSUE_H_
#include <apfSimFEA.h>
namespace bio
{
  class LinearTissue : public amsi::apfSimFEA
  {
    protected:
    std::map<pGEntity, amsi::ElementalSystem*> constitutives;
    public:
    LinearTissue(pGModel imdl, pParMesh imsh, pACase ipd, pACase iss, MPI_Comm cm);
    ~LinearTissue();
    virtual void UpdateDOFs(const double * sol) override;
    virtual void Assemble(amsi::LAS * las) override;
    apf::Field * getField(){ return apf_primary_field; }
  };
}
#endif
