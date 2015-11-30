#include "apfUtil.h"
#include <apfNumbering.h>
#include <cassert>
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
  apf::Vector3 calcLocalCoord(apf::Mesh * msh,
			      apf::MeshEntity * ent,
			      const apf::Vector3 & glb)
  {
    ma::Affine i = ma::getMap(msh,ent);
    return i * glb;
  }
  void calcLocalCoords(apf::DynamicArray<apf::Vector3> & lcl_crds,
		       apf::Mesh * macro_msh,
		       apf::MeshEntity * macro_ent,
		       const apf::DynamicArray<apf::Vector3> & gbl_crds)
  {
    ma::Affine inv = ma::getMap(macro_msh,macro_ent);
    int sz = gbl_crds.getSize();
    lcl_crds.setSize(sz);
    for(int ii = 0; ii < sz; ii++)
      lcl_crds[ii] = inv * gbl_crds[ii];
  }
  void calcDimMeasures(apf::Mesh * msh, int dim, std::vector<double> & msrs)
  {
    apf::MeshEntity * ent = NULL;
    apf::MeshIterator * it = NULL;
    for(it = msh->begin(dim); ent = msh->iterate(it);)
      msrs.push_back(apf::measure(apf::createMeshElement(msh,ent)));
  }
  void getCoords(apf::Mesh * msh,
		 apf::MeshEntity ** vrts,
		 apf::Vector3 * crds,
		 int nm)
  {
    for(int ii = 0; ii < nm; ii++)
      msh->getPoint(vrts[ii],0,crds[ii]);
  }
  double calcDeformedLength(apf::Field * f, apf::MeshEntity * e)
  {
    apf::Mesh * msh = f->getMesh();
    apf::MeshEntity * vs[2];
    msh->getDownward(e,0,&vs[0]);
    apf::Vector3 crds[2];
    getCoords(msh,&vs[0],&crds[0],2);
    apf::Vector3 du[2];
    apf::getVector(f,vs[0],0,du[0]);
    apf::getVector(f,vs[1],0,du[1]);
    return calcDistance(crds[0]+du[0],crds[1]+crds[1]);
  }
  double calcDistance(const apf::Vector3 & a, const apf::Vector3 & b)
  {
    apf::Vector3 c = b - a;
    double v = c * c;
    return sqrt(v);
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
  void fromArray(apf::DynamicVector & to,
     const double * from,
     int sz)
  {
    to.resize(sz);
    to.zero();
    for(int ii = 0; ii < sz; ii++)
      to[ii] = from[ii];
  }
  void fromArray(apf::DynamicMatrix & to,
     const double * from,
     int nr, int nc)
  {
    to.setSize(nr,nc);
    to.zero();
    for(int ii = 0; ii < nr; ii++)
      for(int jj = 0; jj < nc; jj++)
    to(ii,jj) = from[ii*nc + jj];
  }
}
