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
        apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,me);
        //apf::Element * eic = apf::createElement(apf_mesh->getCoordinateField(),mlm);
        //apf::NewArray<apf::Vector3> ic;
        //apf::getVectorNodes(eic,ic);
        //int nds = apf::countNodes(eic);
        apf::Matrix3x3 eps;
        int ip = apf::countIntPoints(mlm,getOrder(mlm));
        for(int ii = 0; ii < ip; ii++)
        {
          apf::getMatrix(strn,me,ii,eps);
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
            nw_hdrs = hdr;
            ++nw_hdrs;
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
            nw_prms = prms;
            ++nw_prms;
            micro_fo_init_data data;
            //for(int jj = 0; jj < nds; jj++)
            //ic[jj].toArray(&data.init_data[jj*3]);
            data.init_data[0] = eps[0][0];
            data.init_data[1] = eps[1][1];
            data.init_data[2] = eps[2][2];
            data.init_data[3] = eps[1][2];
            data.init_data[4] = eps[0][2];
            data.init_data[5] = eps[0][1];
            *nw_data++ = data;
            //++nw_data;
          }
          apf::destroyMeshElement(mlm);
        }
      }
    }
    apf_mesh->end(it);
  }
  template <typename O>
    void MultiscaleTissue::serializeRVEData(O o)
  {
    for(std::list<apf::MeshEntity*>::iterator me = rve_ents.begin(); me != rve_ents.end(); ++me)
    {
      int rve_cnt = countRVEsOn(*me);
      if(rve_cnt > 0)
      {
        apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,*me);
        // apf::Element * e_disp_delta = apf::createElement(delta_u,mlm);
        // apf::Element * e_disp = apf::createElement(apf_primary_field,mlm);
        // apf::Element * e_init_coords = apf::createElement(apf_mesh->getCoordinateField(),mlm);
        rslt_mp[*me].resize(rve_cnt);
        /*
        // disp_deltas holds the solutions of the Newton-Raphson Procedure for the macroscale solve.
        apf::NewArray<apf::Vector3> disp_deltas;
        apf::getVectorNodes(e_disp_delta,disp_deltas);
        // total_disp holds the accumulation of disp_deltas (total amount displaced from intial configuration).
        apf::NewArray<apf::Vector3> total_disp;
        apf::getVectorNodes(e_disp,total_disp);
        apf::NewArray<apf::Vector3> init_coords;
        apf::getVectorNodes(e_init_coords,init_coords);
        int num_nodes = apf::countNodes(e_disp);
        */
        apf::Matrix3x3 eps;
        for(int ii = 0; ii < rve_cnt; ii++)
        {
          apf::getMatrix(strn,*me,ii,eps);
          micro_fo_data fo_data; // fo_data stores macroscale displacement data that is sent to microscale.
          /*
          for(int node = 0; node < num_nodes; node++)
          {
            apf::Vector3 current_coords = init_coords[node] + total_disp[node];
            current_coords.toArray(&fo_data.data[node*3]);
            apf::Vector3 disp_delta = disp_deltas[node];
            disp_delta.toArray(&fo_data.data[node*3 + num_nodes*3]);
          }
          */
          fo_data.data[0] = eps[0][0];
          fo_data.data[1] = eps[1][1];
          fo_data.data[2] = eps[2][2];
          fo_data.data[3] = eps[1][2];
          fo_data.data[4] = eps[0][2];
          fo_data.data[5] = eps[0][1];
          *o++ = fo_data;
        }
        apf::destroyMeshElement(mlm);
      }
    }
  }
}
