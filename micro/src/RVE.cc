#include "RVE.h"
#include "FiberNetwork.h"
#include "las.h"
#include "lasSparskit.h"
namespace bio
{
  void originCenterUnitCube(std::vector<apf::Vector3> & cbe_crnrs)
  {
    // ordering of these loops is important to order the nodes correctly from 0-7
    for(int y = -1; y <= 1; y += 2)
      for(int z = 1; z >= -1; z -= 2)
       for(int x = -1; x <= 1; x += 2)
	 cbe_crnrs.push_back(apf::Vector3(0.5*x,0.5*y,0.5*z));
  }
  void originCenterUnitSquare(std::vector<apf::Vector3> & cbe_crnrs)
  {
    cbe_crnrs.push_back(apf::Vector3(-0.5,-0.5, 0.0));
    cbe_crnrs.push_back(apf::Vector3( 0.5,-0.5, 0.0));
    cbe_crnrs.push_back(apf::Vector3( 0.5, 0.5, 0.0));
    cbe_crnrs.push_back(apf::Vector3(-0.5, 0.5, 0.0));
  }
  RVE::RVE(int d)
  : dim(d)
  , hd(0.5)
  , cbe(NULL)
  , cbe_u_e(NULL)
  , cbe_u(NULL)
  , cbe_dof(NULL)
  {
    assert(d == 2 | d == 3);
    std::vector<apf::Vector3> cbe_crnrs;
    if(dim == 3)
      originCenterUnitCube(cbe_crnrs);
    else if(dim == 2)
      originCenterUnitSquare(cbe_crnrs);
    cbe = makeSingleEntityMesh(apf::Mesh::HEX,&cbe_crnrs[0]);
    apf::MeshIterator * it = cbe->begin(3);
    apf::MeshEntity * cbe_e = cbe->iterate(it);
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
  // technically the vector should only return the nodes, since if an edge
  //  had two mid-edge nodes one of them could theoretically by on the boundary
  //  and the other not, but we will never run into this
  //  also this is basically a FieldOp...
  void calcBoundaryNodes(const RVE * rve,
			 FiberNetwork * fn,
			 std::vector<apf::MeshEntity*> & bnds)
  {
    apf::Mesh * fn_msh = fn->getNetworkMesh();
    int dim = fn->getDim();
    for(int d = 0; d < dim; d++)
    {
      apf::MeshEntity * me = NULL;
      apf::MeshIterator * it = NULL;
      for(it = fn_msh->begin(d); me = fn_msh->iterate(it);)
      {
	int nds = fn->getDisplacementField()->countNodesOn(me);
	for(int nd = 0; nd < nds; nd++)
	{
	  apf::Vector3 crd;
	  fn_msh->getPoint(me,nd,crd);
	  if(rve->onBoundary(crd))
	    bnds.push_back(me);
	}
      }
      fn_msh->end(it);
    }
  }
  void applyRVEForceBC(las::skVec * f,
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