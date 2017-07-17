#include "bioUtil.h"
#include <maAffine.h>
//#include <maAffine.h>
#include <gmi.h>
#include <gmi_null.h>
#include <apfMDS.h>
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
  apf::Matrix3x3 eye()
  {
    apf::Matrix3x3 rslt;
    for(int ii = 0; ii < 3; ii++)
      rslt[ii][ii] = 1.0;
    return rslt;
  }
  apf::Matrix3x3 ones()
  {
    apf::Matrix3x3 rslt;
    for(int ii = 0; ii < 3; ii++)
      for(int jj = 0; jj < 3; jj++)
        rslt[ii][jj] = 1.0;
    return rslt;
  }
  void mapGlobalToLocal(apf::Mesh * msh,
                        apf::MeshEntity * ent,
                        apf::Vector3 const & global,
                        apf::Vector3 & local)
  {
    ma::Affine i = ma::getMap(msh,ent);
    local = i * global;
  }
}
