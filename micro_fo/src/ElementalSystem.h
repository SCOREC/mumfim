#ifndef BIO_ELEMENTAL_SYSTEM_H_
#define BIO_ELEMENTAL_SYSTEM_H_

#include <apf.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>

namespace bio
{
  /**
   * A class to model the elemental system produced by analysis of an individual
   *  element of a mesh. Only models the elemental system, not any relationship between
   *  the elemental system and the global linear system, so dof numbering is not
   *  acknowledged or handled in any way at this level.
   */
  class ElementalSystem
  {
  protected:
    apf::DynamicVector fe;
    apf::DynamicMatrix ke;
  public:
    ElementalSystem()
      : fe()
      , ke()
    { }
    void resize(int nedofs)
    {
      fe.setSize(nedofs);
      ke.setSize(nedofs,nedofs);
    }
    void zero()
    {
      fe.zero();
      ke.zero();
    }
    void addToVector(int idx, double v)
    {
      fe(idx) += v;
    }
    void addToMatrix(int r, int c, double v)
    {
      ke(r,c) += v;
    }
    int nedofs() const
    {
      return fe.getSize();
    }
    const apf::DynamicVector & getfe() const {return fe;}
    const apf::DynamicMatrix & getKe() const {return ke;}
  };


  template <typename A, typename B, typename C>
    void assemble(A from, B indices, C to, int cnt)
  {
    for(int ii = 0; ii < cnt; ii++)
      to[indices[ii]] = from[ii];
  }
}

#endif
