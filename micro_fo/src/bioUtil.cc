#include "bioUtil.h"
#include <maAffine.h>
#include <maMap.h>
namespace bio
{
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
    int tp = msh->getType(ent);
    if(tp == apf::Mesh::EDGE || tp == apf::Mesh::TRIANGLE || tp == apf::Mesh::TET)
    {
      ma::Affine i = ma::getMap(msh,ent);
      local = i * global;
    }
    else
      for(int ii = 0; ii < msh->getDimension(); ++ii)
        local[ii] = 0.0;
    return;
  }
}
