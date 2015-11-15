#ifndef BIO_ELEMENTAL_SYSTEM_H_
#define BIO_ELEMENTAL_SYSTEM_H_

#include <apf.h>

namespace bio
{
  class ElementalSystem
  {
    virtual const apf::DynamicMatrix & getKe() = 0;
    virtual const apf::DynamicVector & getfe() = 0;
  };

  class FacadeElementalSystem : public ElementalSystem
  {
  protected:
    const apf::DynamicMatrix & ke;
    const apf::DynamicVector & fe;
  public:
    FacadeElementalSystem(apf::DynamicMatrix & k, apf::DynamicVector & f)
      : ke(k)
      , fe(f)
    { }
    const apf::DynamicMatrix & getKe() { return ke; }
    const apf::DynamicVector & getFe() { return fe; }
  };

  class SimpleElementalSystem : public ElementalSystem
  {
  protected:
    apf::DynamicMatrix ke;
    apf::DynamicVector fe;
  public:
    
  };
}

#endif
