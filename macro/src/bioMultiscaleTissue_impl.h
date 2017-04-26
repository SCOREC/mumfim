#include "amsiDeformation.h"
namespace bio
{
  template <typename O>
    void MultiscaleTissue::updateRVEDeletion(O o, bool all)
  {
    int iid = 0;
    for(std::list<apf::MeshEntity*>::iterator me = rve_ents.begin();
        me != rve_ents.end(); me++)
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
  template <typename O1, typename O2, typename O3, typename O4>
    void MultiscaleTissue::updateRVEAddition(O1 nw_ents,
                                             O2 nw_hdrs,
                                             O3 nw_prms,
                                             O4 nw_data,
                                             bool all)
  {
    apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
    for(apf::MeshEntity * me = NULL; (me = apf_mesh->iterate(it));)
    {
      // parse simmetrix RVE specification
      pEntity snt = reinterpret_cast<pEntity>(me);
      pGEntity gsnt = EN_whatIn(snt);
      pAttribute mm = GEN_attrib(gsnt,"material model");
      pAttribute sm = Attribute_childByType(mm,"multiscale model");
      if (sm)
      {
        pAttributeTensor0 fbr_rd = (pAttributeTensor0)Attribute_childByType(sm,"radius");
        pAttributeTensor0 vl_frc = (pAttributeTensor0)Attribute_childByType(sm,"volume fraction");
        pAttribute prms = Attribute_childByType(sm,"force reaction");
        pAttributeTensor0 yngs = (pAttributeTensor0)Attribute_childByType(prms,"youngs modulus");
        pAttributeTensor0 nnlr = (pAttributeTensor0)Attribute_childByType(prms,"nonlinearity parameter");
        pAttributeTensor0 lntr = (pAttributeTensor0)Attribute_childByType(prms,"linear transition");
        int fbr_rctn = nnlr ? 1 : 0;
        pAttribute ornt = Attribute_childByType(sm,"fiber orientation");
        bool orntd = ornt != NULL;
        pAttributeTensor1 axs = NULL;
        pAttributeTensor0 algn = NULL;
        if(orntd)
        {
          axs = (pAttributeTensor1)Attribute_childByType(ornt,"axis");
          algn = (pAttributeTensor0)Attribute_childByType(ornt,"alignment");
        }
        // create new data for microscale
        apf::MeshElement * ml = apf::createMeshElement(apf_mesh,me);
        apf::Element * e  = apf::createElement(apf_primary_field,ml);
        apf::Element * ce = apf::createElement(apf_mesh->getCoordinateField(),ml);
        apf::Matrix3x3 F;
        apf::NewArray<apf::Vector3> Ni;
        apf::getVectorNodes(ce,Ni);
        int nds = apf::countNodes(ce);
        int ip = apf::countIntPoints(ml,getOrder(ml));
        for(int ii = 0; ii < ip; ii++)
        {
          apf::Vector3 p;
          apf::getIntPoint(ml,1,ip,p);
          amsi::deformationGradient(e,p,F);
          int crt = apf::getScalar(crt_rve,me,ii);
          int prv = apf::getScalar(prv_rve,me,ii);
          if((crt == FIBER_ONLY && prv != FIBER_ONLY) || all)
          {
            nw_ents = me;
            ++nw_ents;
            micro_fo_header hdr;
            hdr.data[RVE_TYPE]       = getRVEType(reinterpret_cast<apf::ModelEntity*>(gsnt));
            hdr.data[ELEMENT_TYPE]   = apf_mesh->getType(me);
            hdr.data[GAUSS_ID]       = ii;
            hdr.data[FIBER_REACTION] = fbr_rctn;
            hdr.data[IS_ORIENTED]    = orntd;
            *nw_hdrs++ = hdr;
            micro_fo_params prms;
            prms.data[FIBER_RADIUS]    = AttributeTensor0_value(fbr_rd);
            prms.data[VOLUME_FRACTION] = AttributeTensor0_value(vl_frc);
            prms.data[YOUNGS_MODULUS]  = AttributeTensor0_value(yngs);
            prms.data[NONLINEAR_PARAM] = nnlr ? AttributeTensor0_value(nnlr) : 0.0;
            prms.data[LINEAR_TRANSITION] = lntr ? AttributeTensor0_value(lntr) : 0.0;
            prms.data[ORIENTATION_AXIS_X] = orntd ? AttributeTensor1_value(axs,0) : 0.0;
            prms.data[ORIENTATION_AXIS_Y] = orntd ? AttributeTensor1_value(axs,1) : 0.0;
            prms.data[ORIENTATION_AXIS_Z] = orntd ? AttributeTensor1_value(axs,2) : 0.0;
            prms.data[ORIENTATION_ALIGN]  = orntd ? AttributeTensor0_value(algn) : 0.0;
            *nw_prms++ = prms;
            micro_fo_init_data data;
            for(int jj = 0; jj < nds; ++jj)
              Ni[jj].toArray(&data.init_data[jj*3]);
            *nw_data++ = data;
          }
          apf::destroyElement(e);
          apf::destroyMeshElement(ml);
        }
      }
    }
    apf_mesh->end(it);
  }
  // TODO switchover to using deformation gradient means that we need to know which integration point in the element we are talking about when
  //      we serialize the RVE data, previously micro already knew this from initialization so it wasn't important at macro
  template <typename O>
    void MultiscaleTissue::serializeRVEData(O o)
  {
    for(auto me = rve_ents.begin(); me != rve_ents.end(); ++me)
    {
      int rve_cnt = countRVEsOn(*me);
      if(rve_cnt > 0)
      {
        apf::MeshElement * ml = apf::createMeshElement(apf_mesh,*me);
//        apf::Element * e = apf::createElement(apf_primary_field,ml);
	apf::Element * e = apf::createElement(delta_u,ml);
        apf::Matrix3x3 F;
        //int ip = apf::countIntPoints(ml,getOrder(ml));
        rslt_mp[*me].resize(rve_cnt);
        for(int ii = 0; ii < rve_cnt; ii++)
        {
          apf::Vector3 p;
          apf::getIntPoint(ml,1,ii,p);
          amsi::deformationGradient(e,p,F);
          micro_fo_data data;
          for(int ii = 0; ii < 3; ++ii)
            for(int jj = 0; jj < 3; ++jj)
              data.data[ii*3 + jj] = F[ii][jj];
          *o++ = data;
        }
        apf::destroyElement(e);
        apf::destroyMeshElement(ml);
      }
    }
  }
}
