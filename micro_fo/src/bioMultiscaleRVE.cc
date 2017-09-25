#include "bioMultiscaleRVE.h"
#include "bioUtil.h"
#include <apfMeshUtil.h>
#include <cmath> // M_PI
#include <numeric>
#include <cassert>
namespace bio
{
  double calcRVEDimensionality(FiberNetwork * fn,
                               double fbr_area,
                               double fbr_vl_frc)
  {
    std::vector<double> lngths;
    apf::Mesh * fn_msh = fn->getNetworkMesh();
    calcDimMeasures(fn_msh,1,std::back_inserter(lngths));
    double ttl_fbr_lngth = std::accumulate(lngths.begin(),lngths.end(),0.0);
    return sqrt(ttl_fbr_lngth * fbr_area / fbr_vl_frc);
  }
  MultiscaleRVE::MultiscaleRVE(RVE * rve,
                               FiberNetwork * fn,
                               micro_fo_header & hdr,
                               micro_fo_params & prm,
                               micro_fo_init_data & dat)
    : gss_id(hdr.data[GAUSS_ID])
    , dim(rve->getDim())
    , lcl_gss()
    , macro(NULL)
    , macro_ent(NULL)
    , macro_melmnt(NULL)
    , macro_elmnt(NULL)
    , nnd(0)
    , fbr_area(prm.data[FIBER_RADIUS]*prm.data[FIBER_RADIUS] * M_PI)
    , fbr_vl_frc(prm.data[VOLUME_FRACTION])
    , rve_dim(1.0) //assumption
  {
    (void)fbr_area;
    (void)fbr_vl_frc;
    assert(dim == fn->getNetworkMesh()->getDimension());
    apf::getGaussPoint(hdr.data[ELEMENT_TYPE],
                       hdr.data[FIELD_ORDER],
                       gss_id,
                       lcl_gss);
    //rve_dim = calcRVEDimenstionality(fn,fbr_area,fbr_vl_frc);
    rve_dim = 1.0;
    std::vector<apf::Vector3> fe_nds;
    int nv = apf::Mesh::adjacentCount[hdr.data[ELEMENT_TYPE]][0];
    for(int ii = 0; ii < nv; ++ii)
      fe_nds.push_back(apf::Vector3(dat.init_data[ii*3],dat.init_data[ii*3+1],dat.init_data[ii*3+2]));
    macro = amsi::makeSingleEntityMesh(apf::Mesh::Type(hdr.data[ELEMENT_TYPE]),&fe_nds[0]);
    apf::MeshIterator * it = macro->begin(3);
    macro_ent = macro->iterate(it);
    macro->end(it);
    macro_u = apf::createLagrangeField(macro,"u",apf::VECTOR,hdr.data[FIELD_ORDER]);
    apf::zeroField(macro_u);
    macro_melmnt = apf::createMeshElement(macro,macro_ent);
    macro_elmnt = apf::createElement(macro_u,macro_ent);
    apf::FieldShape * shp = apf::getShape(macro_u);
    nnd = shp->countNodesOn(hdr.data[ELEMENT_TYPE]);
  }
  MultiscaleRVE::~MultiscaleRVE()
  {
    apf::destroyElement(macro_elmnt);
    apf::destroyMeshElement(macro_melmnt);
  }
  void MultiscaleRVE::dCidFE(apf::DynamicMatrix & dRVEdFE,
                             const int ii,
                             const apf::Vector3 & ci,
                             const double rve_dim)
  {
    apf::NewArray<double> N;
    apf::getShapeValues(macro_elmnt,ci,N);
    for(int jj = 0; jj < nnd; ++jj)
      for(int kk = 0; kk < dim; ++kk)
        dRVEdFE(ii*kk,jj*kk) = N[jj] / rve_dim;
  }
  void MultiscaleRVE::calcdRVEdFE(apf::DynamicMatrix & dRVEdFE,
                                  const RVE * rve)
  {
    int num_rve_nds = rve->numNodes();
    apf::Vector3 gbl_gss;
    apf::mapLocalToGlobal(macro_melmnt,lcl_gss,gbl_gss);
    apf::DynamicArray<apf::Vector3> gbl_rve_crds(num_rve_nds);
    calcGlobalRVECoords(gbl_rve_crds,rve_dim,gbl_gss);
    std::vector<apf::Vector3> lcl_rve_crds(num_rve_nds);
    mapGlobalsToLocals(macro,
                       macro_ent,
                       gbl_rve_crds.begin(),
                       gbl_rve_crds.end(),
                       std::back_inserter(lcl_rve_crds));
    dRVEdFE.zero();
    int sz = gbl_rve_crds.size();
    for(int ii = 0; ii < sz; ++ii)
      dCidFE(dRVEdFE,ii,lcl_rve_crds[ii]-lcl_gss,dim);
  }
}
