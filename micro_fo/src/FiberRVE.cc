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

    std::vector<apf::Vector3> cbe_crnrs;
    // ordering of these loops is important to order the nodes correctly from 0-7
    for(int y = -1; y <= 1; y += 2)
      for(int z = -1; z <= 1; z += 2)
	for(int x = -1; x <= 1; x += 2)
	  cbe_crnrs.push_back(apf::Vector3(0.5*x,0.5*y,0.5*z));
    cbe = makeSingleEntityMesh(apf::Mesh::HEX,&cbe_crnrs[0]);
  }

  void FiberRVE::forwardCubeDisp()
  {
    apf::MeshIterator * it = cbe->begin(dim);
    apf::MeshEntity * cbe_me = cbe->iterate(it);
    cbe->end(it);
    apf::Element * cbe_e = apf::createElement(cbe_u,apf::createMeshElement(cbe,cbe_me));
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
}
