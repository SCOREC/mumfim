#include "bioMultiscaleRVE.h"
#include "bioUtil.h"
#include <apfMeshUtil.h>
#include <apfMDS.h>
#include <apfConvert.h>
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
  MultiscaleRVE::MultiscaleRVE(const MultiscaleRVE& mrve)
  {
    gss_id = mrve.gss_id;
    dim = mrve.dim;
    lcl_gss = mrve.lcl_gss;
    macro = amsi::makeNullMdlEmptyMesh();
    apf::convert(mrve.macro, static_cast<apf::Mesh2*>(macro));
    macro_u = macro->findField(apf::getName(mrve.macro_u));
    apf::MeshIterator* it = macro->begin(3);
    macro_ent = macro->iterate(it);
    macro->end(it);
    macro_melmnt = apf::createMeshElement(macro, macro_ent);
    macro_elmnt = apf::createElement(macro_u, macro_ent);
    nnd = mrve.nnd;
    fbr_area = mrve.fbr_area;
    rve_dim = mrve.rve_dim;
    scale_conversion = mrve.scale_conversion;
    fbr_vl_frc = mrve.fbr_vl_frc;
  }
  MultiscaleRVE::MultiscaleRVE(RVE* r, FiberNetwork* fn, micro_fo_header& hdr,
                               micro_fo_params& prm, micro_fo_init_data& dat)
      : rve(r)
      , gss_id(hdr.data[GAUSS_ID])
      , dim(rve->getDim())
      , lcl_gss()
      , macro(NULL)
      , macro_ent(NULL)
      , macro_melmnt(NULL)
      , macro_elmnt(NULL)
      , nnd(0)
      , fbr_area(prm.data[FIBER_RADIUS] * prm.data[FIBER_RADIUS] * M_PI)
      , fbr_vl_frc(prm.data[VOLUME_FRACTION])
      , rve_dim(1.0)
      , scale_conversion()
  {
    (void)fbr_area;
    (void)fbr_vl_frc;
    assert(dim == fn->getNetworkMesh()->getDimension());
    apf::getGaussPoint(hdr.data[ELEMENT_TYPE],
                       hdr.data[FIELD_ORDER],
                       gss_id,
                       lcl_gss);
    rve_dim = calcRVEDimensionality(fn,fbr_area,fbr_vl_frc);
    scale_conversion = 1.0 / (rve_dim * rve_dim);
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
    nnd = apf::countNodes(macro_elmnt);
  }
  MultiscaleRVE::~MultiscaleRVE()
  {
    apf::destroyField(macro_u);
    apf::destroyElement(macro_elmnt);
    apf::destroyMeshElement(macro_melmnt);
    macro->destroyNative();
    apf::destroyMesh(macro);
  }
  void MultiscaleRVE::calcdRVEdFE(apf::DynamicMatrix & dRVEdFE)
  {
    int num_rve_nds = rve->numNodes();
    apf::NewArray<int> rw_ids;
    apf::getElementNumbers(rve->getNumbering(),rve->getMeshEnt(),rw_ids);
    apf::Vector3 gbl_gss;
    // calculated from the reference element, 
    //  not the displaced element
    apf::mapLocalToGlobal(macro_melmnt,lcl_gss,gbl_gss);
    apf::DynamicArray<apf::Vector3> gbl_rve_crds(num_rve_nds);
    calcGlobalRVECoords(rve,gbl_rve_crds,rve_dim,gbl_gss);
    std::vector<apf::Vector3> lcl_rve_crds;
    mapGlobalsToLocals(macro,
                       macro_ent,
                       gbl_rve_crds.begin(),
                       gbl_rve_crds.end(),
                       std::back_inserter(lcl_rve_crds));
    int sz = gbl_rve_crds.size();
    dRVEdFE.setSize(sz*dim,nnd*dim);
    dRVEdFE.zero();
    for(int ii = 0; ii < sz; ++ii)
    {
      apf::NewArray<double> gss_shps;
      apf::getShapeValues(macro_elmnt,lcl_gss,gss_shps);
      apf::NewArray<double> rvec_shps;
      apf::getShapeValues(macro_elmnt,lcl_rve_crds[ii],rvec_shps);
      for(int jj = 0; jj < nnd; jj++)
      {
        double diff  = rvec_shps[jj] - gss_shps[jj];
        diff /= rve_dim;
        for(int kk = 0; kk < dim; ++kk)
          dRVEdFE(rw_ids[ii*dim+kk],jj*dim+kk) = diff;
      }
    }
  }
}
