#include "RVE.h"
#include "MacroCoupling.h"

#include <apf.h>
#include <apfDynamicMatrix.h>

namespace bio
{
  MacroInfo::MacroInfo()
    : gss_id()
    , dim()
    , lcl_gss()
    , macro_msh(NULL)
    , macro_ent(NULL)
    , macro_melmnt(NULL)
    , macro_elmnt(NULL)
    , nnd()
    , fbr_area()
    , fbr_vl_frc()
    , rve_dim()
  {}
  
  void MacroInfo::dCidFE(apf::DynamicMatrix & dRVEdFE, const int ii, const apf::Vector3 & ci,
			 double rve_dim)
  {
    apf::NewArray<double> N;
    apf::getShapeValues(macro_elmnt,ci,N);
    for(int jj = 0; jj < nnd; jj++) 
      for(int kk = 0; kk < dim; kk++)
	dRVEdFE(ii*kk,jj*kk) = N[jj] / rve_dim;
  }

  void MacroInfo::calcdRVEdFE(apf::DynamicMatrix & drve_dfe, const RVE * rve)
  {
    int num_rve_nds = rve->numNodes();
      
    apf::Vector3 gbl_gss;
    apf::mapLocalToGlobal(macro_melmnt,lcl_gss,gbl_gss);
    apf::DynamicArray<apf::Vector3> gbl_rve_crds(num_rve_nds);
    calcGlobalRVECoords(gbl_rve_crds,rve_dim,gbl_gss);
    
    apf::DynamicArray<apf::Vector3> lcl_rve_crds(num_rve_nds);
    calcLocalCoords(lcl_rve_crds,macro_msh,macro_ent,gbl_rve_crds);
    
    drve_dfe.zero();
    int sz = lcl_rve_crds.getSize();
    for(int ii = 0; ii < sz; ii++)
      dCidFE(drve_dfe,ii,lcl_rve_crds[ii]-lcl_gss,dim);
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
