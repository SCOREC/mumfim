#include "RVE.h"

namespace bio
{
  RVE::RVE()
    : dim()
    , hd(0.5)
    , cbe(NULL)
    , cbe_u_e(NULL)
    , cbe_u(NULL)
  {
    std::vector<apf::Vector3> cbe_crnrs;
    // ordering of these loops is important to order the nodes correctly from 0-7
    for(int y = -1; y <= 1; y += 2)
      for(int z = 1; z >= -1; z -= 2)
	for(int x = -1; x <= 1; x += 2)
	  cbe_crnrs.push_back(apf::Vector3(hd*x,hd*y,hd*z));
    cbe = makeSingleEntityMesh(apf::Mesh::HEX,&cbe_crnrs[0]);

    //make cbe_u and cbe_e_u
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

  void calcBoundaryNodes(const RVE * rve,
			 const FiberNetwork * fn,
			 std::vector<apf::MeshEntity*> & bnds)
  {
    const apf::Mesh * fn_msh = fn->getNetworkMesh();
    int dim = fn->getDim();
    for(int d = 0; d < dim; d++)
    {
      apf::MeshEntity * me = NULL;
      apf::MeshIterator * it = NULL;
      for(it = fn_msh->begin(d); me = fn_msh->iterate(it);)
      {
	apf::Vector3 crd;
	fn_msh->getPoint(me,0,crd);
	if(rve->onBoundary(crd))
	  bnds.push_back(me);
      }
      fn_msh->end(it);
    }
  }

  void applyRVEForceBC(LinearSystem * ls,
		       RVE * rve,
		       FiberNetwork * fn)
  {
    const apf::Numbering * nm = fn->getNumbering();
    const apf::Field * u = fn->getDisplacementField();

    Vector * vec = ls->getVector();
    LSOps * ops = ls->getOps();
    
    std::vector<apf::MeshEntity*> bnds;
    calcBoundaryNodes(rve,fn,bnds);
    for(std::vector<apf::MeshEntity*>::iterator it = bnds.begin(); it != bnds.end(); it++)
    {
      apf::NewArray<int> dofs;
      apf::getElementNumbers(nm,*it,dofs);
      int nedofs = apf::countNodes(apf::createElement(u,*it));
      apf::DynamicVector zeros(nedof);
      zeros.zero();
      ops->setVector(vec,zeros,dofs,nedofs);
    }
  }
}
