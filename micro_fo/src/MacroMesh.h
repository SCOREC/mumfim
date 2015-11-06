#ifndef BIO_MACRO_COUPLING_H
#define BIO_MACRO_COUPLING_H

#include "apfUtil.h"

namespace bio
{
  
  class MacroInfo
  {
  protected:
    int gss_id;
    apf::Mesh * macro_msh;
    apf::MeshEntity * macro_ent;
    int nnd; // num nodes on macro element
  public:
    
  };
  
}

#endif
