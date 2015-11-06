#ifndef BIO_MACRO_COUPLING_H
#define BIO_MACRO_COUPLING_H

#include "apfUtil.h"
#include <apf.h>
#include <apfDynamicMatrix.h>

namespace bio
{
  class MacroInfo
  {
  protected:
    int gss_id;
    apf::Vector3 lcl_gss;
    apf::Mesh * macro_msh;
    apf::MeshEntity * macro_ent;
    apf::MeshElement * macro_melmnt;
    apf::Element * macro_elmnt;
    int nnd; // num nodes on macro element

    void dCidFE(apf::DynamicMatrix&,const int,const apf::Vector3 &, const apf::Vector3 &);
  public:
    void calcdRVEdFE(apf::DynamicMatrix & drve_dfe, const FiberRVE * rve);
  };
  
}

#endif
