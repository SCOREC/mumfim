#include "amsiDeformation.h"
namespace bio
{
  /*
    template <typename O>
    void MultiscaleTissue::updateRVEDeletion(O o, bool all)
    {
    int iid = 0;
    for(auto me = rve_ents.begin(); me != rve_ents.end(); me++)
    {
    apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,*me);
    int num_gauss_pts = apf::countIntPoints(mlm,getOrder(mlm));
    for(int ii = 0 ; ii < num_gauss_pts; ii++)
    {
    int crt = apf::getScalar(crt_rve,*me,ii);
    int prv = apf::getScalar(prv_rve,*me,ii);
    if((crt == NONE && prv != NONE) || all)
    {
    o = iid;
    ++o;
    // assumed that removing one RVE from an entity removed that entity from the list of entities containing RVEs, not necessarily the case when using multiply IP per ent
    //me = rve_ents.erase(me);
    //me--;
    }
    iid++;
    }
    apf::destroyMeshElement(mlm);
    }
    }
  */
  template <typename O1, typename O2, typename O3, typename O4, typename O5>
  void MultiscaleTissue::serializeNewRVEData(O1 new_hdrs, O2 new_prms,
                                             O3 new_data, O4 new_slvr_prms,
                                             O5 new_int_slvr_prms, bool all)
  {
    apf::MeshEntity * rgn = NULL;
    apf::MeshIterator * it = apf_mesh->begin(3);
    while(rgn = apf_mesh->iterate(it))
    {
      apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,rgn);
      int ng = apf::countIntPoints(mlm,getOrder(mlm));
      for(int ip = 0; ip < ng; ++ip)
      {
        MicroscaleType crt = static_cast<MicroscaleType>(apf::getScalar(crt_rve,rgn,ip));
        MicroscaleType prv = static_cast<MicroscaleType>(apf::getScalar(prv_rve,rgn,ip));
        // if the RVE is new
        if((crt == MicroscaleType::FIBER_ONLY && prv != MicroscaleType::FIBER_ONLY) ||
           (crt == MicroscaleType::ISOTROPIC_NEOHOOKEAN && prv != MicroscaleType::ISOTROPIC_NEOHOOKEAN)
            || all) 
        {
          micro_fo_header hdr;
          micro_fo_params prm;
          micro_fo_init_data dat;
          micro_fo_solver slvr;
          micro_fo_int_solver int_slvr;
          getInternalRVEData(rgn,hdr,prm,dat);
          getExternalRVEData(rgn,hdr,prm,slvr,int_slvr);
          hdr.data[GAUSS_ID] = ip;
          *new_hdrs++ = hdr;
          *new_prms++ = prm;
          *new_data++ = dat;
          *new_slvr_prms++ = slvr;
          *new_int_slvr_prms++ = int_slvr;
        }
      }
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(it);
  }
  // TODO switchover to using deformation gradient means that we
  //  oneed to know which integration point in the element we are
  //  talking about when we serialize the RVE data, previously
  //  micro already knew this from initialization so it wasn't
  //  important at macro
  template <typename O>
    void MultiscaleTissue::serializeRVEData(O o)
  {
    apf::MeshEntity * rgn = NULL;
    for(auto * it = apf_mesh->begin(3); (rgn = apf_mesh->iterate(it)); )
    {
      apf::MeshElement * mlm;
      apf::Element * e;
      // on the first time through use the total deformation gradient
      // to capture any deformation from initialization/guess phase of
      // solution which has not been applied to the microscale yet
      if(load_step == 0 && iteration == 0)
      {
        mlm = apf::createMeshElement(apf_mesh,rgn);
        e = apf::createElement(apf_primary_field,mlm);
      }
      else
      {
        mlm = apf::createMeshElement(prev_coords,rgn);
        e = apf::createElement(delta_u,mlm);
      }
      int ng = apf::countIntPoints(mlm,getOrder(mlm));
      for(int ip = 0; ip < ng; ++ip)
      {
        MicroscaleType crt = static_cast<MicroscaleType>(apf::getScalar(crt_rve,rgn,ip));
        if(crt != MicroscaleType::NONE)
        {
          apf::Matrix3x3 F;
          apf::Vector3 p;
          apf::getIntPoint(mlm,1,ip,p);
          amsi::deformationGradient(e,p,F);
          // DEBUG
          if(load_step == 0 && iteration == 0)
          {
            apf::setMatrix(test_inc_dfm, rgn, ip, F);
          }
          else
          {
            apf::Matrix3x3 prv_F;
            apf::getMatrix(test_inc_dfm, rgn, ip, prv_F);
            apf::setMatrix(test_inc_dfm, rgn, ip, F*prv_F);
          }
          // DEBUG
          micro_fo_data data;
          for(int ii = 0; ii < 3; ++ii)
            for(int jj = 0; jj < 3; ++jj)
              data.data[ii*3 + jj] = F[ii][jj];
          *o++ = data;
        }
      }
      apf::destroyElement(e);
      apf::destroyMeshElement(mlm);
    }
  }
}
