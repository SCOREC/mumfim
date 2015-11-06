#include "FiberRVE.h"

#include <apf.h>
#include <maMap.h>

namespace bio
{
  // assumes 3d
  void calcGlobalRVECorners(const apf::Vector3 & coord, double half_rve_dim, apf::Vector3 (&cs)[num_corners])
  {
    static double op[][] = {{-1.0,-1.0, 1.0},
			    { 1.0,-1.0, 1.0},
			    {-1.0,-1.0,-1.0},
			    { 1.0 -1.0,-1.0},
			    {-1.0, 1.0, 1.0},
			    { 1.0, 1.0, 1.0},
			    {-1.0, 1.0,-1.0},
			    { 1.0, 1.0,-1.0}};
    
    for(int ii = 0; ii < num_corners; ii++}
      for(int jj = 0; jj < dim; jj++)
	cs[ii][jj] = coord[jj] + half_rve_dim*op[ii][jj]
  }

  void FiberRVE::dCidFE(apf::DynamicMatrix * dRVEdFE,
			int ii, const apf::Vector & ci, int nsf)
  {
    apf::NewArray<double> N;
    apf::getShapeValues(e,d,N);
    for(int jj = 0; jj < nsf; jj++) 
      for(int kk = 0; kk < dim; kk++)
	dRVEdFE(ii*kk  ,jj*kk  ) = N[jj] / rve_dim;
  }


  // make stateless, pass in RVE
  void FiberRVE::calcdRVEdFE(apf::DynamicMatrix & dRVEdFE)
  {
    dRVEdFE.zero();
    apf::Vector3 g;  // gauss coords
    apf::getGaussPt(apf::getEntityType(me),1,gid,g);
    apf::Vector3 gg; // global gauss coords 
    apf::mapLocalToGlobal(me,g,gg);
    apf::Vector3 cnrs[num_corners];
    calcGlobalRVECorners(gg,half_rve_dim,gbl_cnrs);
    ma::Affine inv = ma::getMap(macro_msh,macro_ent);
    std::transform(&cnrs[0],&cnrs[num_corners],&cnrs[0],std::bind1st(apply,inv));
    int nsf = apf::countNodes(apf::getLagrange(0),macro_msh->getType(macro_ent));
    for(int ii = 0; ii < num_corners; ii++) //rve corners will always be 8...
      dCidFE(dRVEdFE,ii,cnrs[ii]-q,nsf);
  }


  void FiberRVE::genInitGuess()
  {
    // cosider pulling nsf into the clas.s.?
    int nsf = apf::countNodes(apf::getLagrange(0),macro_msh->getType(macro_ent));
    apf::DynamicMartix dRVEdFE(dim*num_corners,dim*nsf);
    calcdRVEdFE(dRVEdFE);
    apf::DynamicVector uFE(dim*nsf);
    apf::DynamicVector uRVE = dRVEdFE * uFE;

    // generate initial guess 
  }

  apf::Vector3 FiberRVE::trilinearInterp(const apf::Vector3 & crd,
					 apf::DynamicMatrix & uRVE)
  {
    return apf::Vector3(0.0,0.0,0.0);
  }

  apf::Vector3 FiberRVE::firstOrderCont(const apf::Vector3 & crd,
					apf::DynamicMatrix & uRVE)
  {
    // all displacements
    //apf::DynamicVector u = dydxr * uRVE;
    retrn apf::Vector3(0.0,0.0,0.0);
  }

  // maybe just build a bi-unit cubic mesh and use the lagrange shape functions to interpolate the displacement from the corners to the nodes?
  void FiberRVE::initialDisplace(apf::MeshEntity * v,
				 const apf::DynamicMatrix & uRVE,
				 apf::Vector3 (*mthd)(const apf::Vector3&, apf::DynamicMatrix&))
  {
    apf::Vector3 c;
    macro_msh->getPoint(v,c);
    apf::Vector3 u = (*mthd)(c,uRVE);
    apf::setVector(micro_u,v,0,u);
  }
}
