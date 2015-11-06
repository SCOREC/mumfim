#ifndef BIO_APF_UTIL_H_
#define BIO_APF_UTIL_H_

#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_null.h>

namespace bio
{

  apf::Mesh2 * makeNullMdlEmptyMesh()
  {
    gmi_register_null();
    gmi_model * mdl = gmi_load(".null");
    return apf::makeEmptyMdsMesh(mdl,3,false);
  };

  apf::Mesh2 * makeSingleEntityMesh(apf::Mesh::Type t, const apf::Vector3 * vs)
  {
    apf::Mesh2 * msh = makeNullMdlEmptyMesh();
    apf::buildOneElement(msh,NULL,t,vs);
    msh->acceptChanges();
    return msh;
  }
}

#endif
