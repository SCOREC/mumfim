#include "bioUtil.h"
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
  class Affine
  {
  public:
    apf::Vector3 operator*(apf::Vector3 const& x) const
    {
      return A * x + b;
    }
    apf::Matrix3x3 A;
    apf::Vector3 b;
  };
  Affine getTetMap(apf::Mesh * msh, apf::MeshEntity * tet)
  {
    apf::Downward v;
    msh->getDownward(tet,0,v);
    apf::Vector3 p[4];
    for(int ii = 0; ii < 4; ++ii)
      msh->getPoint(v[ii],0,p[ii]);
    Affine a;
    a.A[0] = p[1] - p[0];
    a.A[1] = p[2] - p[0];
    a.A[2] = p[3] - p[0];
    a.A = apf::transpose(a.A);
    a.A = apf::invert(a.A);
    a.b = p[0];
    a.b = (a.A * a.b) * -1;
    return a;
  }
  // only works for linear tets.
  void mapGlobalToLocal(apf::Mesh * msh,
                        apf::MeshEntity * ent,
                        apf::Vector3 const & global,
                        apf::Vector3 & local)
  {
    (void)global;
    int tp = msh->getType(ent);
    if(tp == apf::Mesh::TET)
    {
      Affine i = getTetMap(msh,ent);
      local = i * global;
    }
    else
      for(int ii = 0; ii < msh->getDimension(); ++ii)
        local[ii] = 0.0;
    return;
  }
}
