#ifndef BIO_MACRO_MESH_H
#define BIO_MACRO_MESH_H

#include "apfUtil.h"

namespace bio
{
  apf::Mesh2 * makeSingleEntityMesh(apf::Mesh::Type t, const apf::Vector3 * vs)
  {
    apf::Mesh2 * msh = makeNullMdlEmptyMesh();
    msh->buildOneElement(msh,NULL,t,vs);
    msh->acceptChanges();
    return msh;
  }
}

#endif
