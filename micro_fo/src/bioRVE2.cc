#include "RVE.h"
#include "FiberNetwork.h"
#include "las.h"
#include "lasSparskit.h"
namespace bio
{
  template <typename O>
  void originCenterCube(O cbe_crnrs, double crd)
  {
    // ordering of these loops is important to order the nodes correctly from 0-7
    for(int y = -1; y <= 1; y += 2)
      for(int z = 1; z >= -1; z -= 2)
       for(int x = -1; x <= 1; x += 2)
	 *cbe_crnrs++ = apf::Vector3(crd*x,crd*y,crd*z);
  }
  template <typename O>
  void originCenterSquare(O cbe_crnrs, double crd)
  {
    *cbe_crnrs++ = apf::Vector3(-crd,-crd, 0.0);
    *cbe_crnrs++ = apf::Vector3( crd,-crd, 0.0);
    *cbe_crnrs++ = apf::Vector3( crd, crd, 0.0);
    *cbe_crnrs++ = apf::Vector3(-cdr, crd, 0.0);
  }
  RVE::RVE(int d)
  : dim(d)
  , crd(0.5)
  , cbe(NULL)
  , cbe_e(NULL)
  , cbe_u_e(NULL)
  , cbe_u(NULL)
  , cbe_dof(NULL)
  {
    assert(d == 2 | d == 3);
    std::vector<apf::Vector3> cbe_crnrs;
    if(dim == 3)
      originCenterUnitCube(cbe_crnrs,hd);
    else if(dim == 2)
      originCenterUnitSquare(cbe_crnrs,hd);
    cbe = makeSingleEntityMesh(apf::Mesh::HEX,&cbe_crnrs[0]);
    apf::MeshIterator * it = cbe->begin(3);
    cbe_e = cbe->iterate(it);
    cbe->end(it);
    cbe_u = apf::createLagrangeField(cbe,"u",apf::VECTOR,1);
    cbe_dof = apf::createNumbering(cbe_u);
    cbe_u_e = apf::createElement(cbe_u,cbe_e);
  }
  void forwardRVEDisplacement(RVE * rve, FiberNetwork * fn)
  {
    apf::Element * cbe_e = rve->getElement();
    apf::Field * fn_u = fn->getDisplacementField();
    InterpOnto forward_u(cbe_e,fn_u);
    forward_u.apply(fn_u);
  }
  void calcGlobalRVECoords(apf::DynamicArray<apf::Vector3> & rve_crds,
			   double rve_dim,
			   const apf::Vector3 & gbl_gss)
  {
    static double op[8][3] = {{-1.0,-1.0, 1.0},
			      { 1.0,-1.0, 1.0},
			      {-1.0,-1.0,-1.0},
			      { 1.0 -1.0,-1.0},
			      {-1.0, 1.0, 1.0},
			      { 1.0, 1.0, 1.0},
			      {-1.0, 1.0,-1.0},
			      { 1.0, 1.0,-1.0}};
    int sz = rve_crds.getSize();
    int d = sz == 8 ? 3 : 2;
    int o = d == 3 ? 0 : 2;
    for(int ii = 0; ii < sz; ii++)
      for(int jj = 0; jj < d; jj++)
	rve_crds[ii][jj] = gbl_gss[jj] + rve_dim*op[ii+o][jj];
  }
  void applyRVEForceBC(las::Vec * f,
		       RVE * rve,
		       FiberNetwork * fn)
  {
    apf::Numbering * nm = fn->getNumbering();
    apf::Field * du = fn->getDisplacementField();
    std::vector<apf::MeshEntity*> bnds;
    calcBoundaryNodes(rve,fn,bnds);
    for(std::vector<apf::MeshEntity*>::iterator it = bnds.begin(); it != bnds.end(); it++)
    {
      apf::NewArray<int> dofs;
      int nedofs = apf::getElementNumbers(nm,*it,dofs);
      apf::DynamicVector zeros(nedofs);
      zeros.zero();
      las::setVecValues(f,zeros,dofs,nedofs,false);
    }
  }
  void displaceRVE(RVE * rve,const apf::DynamicVector & du)
  {
    ApplySolution(rve->getNumbering(),&du[0],0,true).apply(rve->getField());
  }
}
