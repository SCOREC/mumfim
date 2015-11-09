#include "FiberRVE.h"
#include "InitGuess.h"
#include "apfUtil.h"

#include <apfShape.h>
#include <apfMesh.h>

#include <cmath>
#include <numeric>

namespace bio
{
  FiberRVE::FiberRVE(apf::Mesh * f) :
    fbr_area(0.0),
    fbr_vl_frc(0.003),
    rve_dim(0.0),
    cbe(NULL),
    cbe_u(NULL),
    fn(f),
    fn_u(NULL),
    dim(f->getDimension())
  {
    double fbr_rds = 3.49911271e-8;
    
    std::vector<double> fbr_lngths;
    calcEdgeLengths(fn,fbr_lngths);
    double ttl_fbr_lngth = std::accumulate(fbr_lngths.begin(),fbr_lngths.end(),0.0);
    fbr_area = fbr_rds * fbr_rds * M_PI;
    rve_dim = sqrt(ttl_fbr_lngth * fbr_area /  fbr_vl_frc);
  }
  
  void FiberRVE::interpCornerDisps()
  {
    apf::Vector3 cu[8];
    // cu - corner displacements
    double cfs[8] = {};
    apf::Vector3 dms(rve_dim,rve_dim,rve_dim);
    apf::Vector3 trn(-0.5,-0.5,-0.5);
    for(int d = 0; d < dim; d++)
    {
      apf::MeshEntity * me = NULL;
      apf::MeshIterator * it = NULL;
      for(it = fn->begin(d); me = fn->iterate(it); )
      {
	apf::FieldShape * fs = apf::getShape(fn_u);
	int nds = fs->countNodesOn(fn->getType(me));
	for(int nd = 0; nd < nds; nd++)
	{
	  tricubicInterpCoefs(cfs,dms,trn,cfs);
	  apf::Vector3 u;
	  for(int ii = 0; ii < 8; ii++)
	    u += cu[ii] * cfs[ii];
	}
      }
      fn->end(it);
    }
  }

  void FiberRVE::forwardCubeDisp()
  {
    apf::MeshIterator * it = cbe->begin(dim);
    apf::MeshEntity * cbe_me = cbe->iterate(it);
    cbe->end(it);
    apf::Element * cbe_e = apf::createElement(cbe_u,apf::createMeshElement(cbe,cbe_me));
    
    for(int d = 0; d < dim; d++)
    {
      NODE_ITER(d,fn,fn_u,
		apf::Vector3 p;
		apf::Vector3 cbe_u;
		fn->getPoint(me,nd,p);
		apf::getVector(cbe_e,p,cbe_u);
		apf::setVector(fn_u,me,nd,cbe_u);
	);
    }
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
}
