#ifndef BIO_MASS_INTEGRATOR_H_
#define BIO_MASS_INTEGRATOR_H_
namespace bio
{
  class MassIntegrator : public apf::integrator
  {
    protected:
      apf::Field * density;
      apf::Numbering * nm;
      apf::Element* elmt;

  };
}

#endif
