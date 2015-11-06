#include "FiberRVE.h"
#include "MacroCoupling.h"

#include <apf.h>
#include <apfDynamicMatrix.h>
#include <maMap.h>

namespace bio
{
  // move to FiberRVE.cc
  void calcGlobalRVECoords(apf::DynamicArray<apf::Vector3> & rve_crds,
			   const apf::Vector3 & rve_dim,
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
	rve_crds[ii][jj] = gbl_gss[jj] + rve_dim[jj]*op[ii+o][jj];
  }

  // move to apf util i guess
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
  
  void MacroInfo::dCidFE(apf::DynamicMatrix & dRVEdFE, const int ii, const apf::Vector3 & ci,
			 const apf::Vector3 & rve_dim)
  {
    apf::NewArray<double> N;
    apf::getShapeValues(macro_elmnt,ci,N);
    for(int jj = 0; jj < nnd; jj++) 
      for(int kk = 0; kk < dim; kk++)
	dRVEdFE(ii*kk,jj*kk) = N[jj] / rve_dim[kk];
  }

  void MacroInfo::calcdRVEdFE(apf::DynamicMatrix & drve_dfe, const FiberRVE * rve)
  {
    int num_rve_nds = rve->numCornerNodes();
    apf::Vector3 rve_dims = rve->getRVEDimensions();
    
    apf::Vector3 gbl_gss;
    apf::mapLocalToGlobal(macro_melmnt,lcl_gss,gbl_gss);
    apf::DynamicArray<apf::Vector3> gbl_rve_crds(num_rve_nds);
    calcGlobalRVECoords(gbl_rve_crds,rve_dims,gbl_gss);
    apf::DynamicArray<apf::Vector3> lcl_rve_crds(num_rve_nds);
    calcLocalCoords(lcl_rve_crds,macro_msh,macro_ent,gbl_rve_crds);
    
    drve_dfe.zero();
    int sz = lcl_rve_crds.getSize();
    for(int ii = 0; ii < sz; ii++)
      dCidFE(drve_dfe,ii,lcl_rve_crds[ii]-lcl_gss,rve_dims);
  }

  /*
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
  */
}
