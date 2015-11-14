#ifndef BIO_APF_UTIL_H_
#define BIO_APF_UTIL_H_

#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_null.h>
#include <maMap.h>

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

  apf::Vector3 calcLocalCoord(apf::Mesh * msh, apf::MeshEntity * ent,const apf::Vector3 & glb)
  {
    ma::Affine i = ma::getMap(msh,ent);
    return i * glb;
  }

  void calcLocalCoords(apf::DynamicArray<apf::Vector3> & lcl_crds,
		       apf::Mesh * macro_msh,
		       apf::MeshEntity * macro_ent,
		       apf::DynamicArray<apf::Vector3> & gbl_crds)
  {
    ma::Affine inv = ma::getMap(macro_msh,macro_ent);
    int sz = gbl_crds.getSize();
    lcl_crds.setSize(sz);
    for(int ii = 0; ii < sz; ii++)
      lcl_crds[ii] = inv * gbl_crds[ii];
  }

  void calcLocalCoord(apf::Mesh * msh,
		      apf::MeshEntity * me,
		      const apf::Vector3 & xyz,
		      apf::Vector3 & xi)
  {
    ma::Affine inv = ma::getMap(msh,me);
    xi =  inv * xyz;
  }

  void calcEdgeLengths(apf::Mesh * msh, std::vector<double> & lngths)
  {
    apf::MeshEntity * edge = NULL;
    apf::MeshIterator * it = NULL;
    for(it = msh->begin(1); edge = msh->iterate(it);)
      lngths.push_back(apf::measure(apf::createMeshElement(msh,edge)));
  }

  void calcEdgeVertDiffs(apf::Mesh * msh, apf::MeshEntity * me, apf::Vector3 & diffs)
  {
    apf::MeshEntity * vs[2];
    msh->getDownward(me,0,&vs[0]);
    apf::Vector3 c[2];
    msh->getPoint(vs[0],0,c[0]);
    msh->getPoint(vs[1],0,c[1]);
    diffs = c[1] - c[0];
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
}

#endif
