#include "bioStiffnessVariation.h"
namespace bio
{
  StiffnessVariation * buildStiffnessVariation(pAcase cs, pANode nd, apf::Field * stf_vrtn)
  {
    StiffnessVariation * stfvar = NULL;

    /* cs = "problem definition" in simmodeler. We are interested in three of teh children.
       - "stiffness source"   : pointer to model region tags.
       - "stiffness sink"     : pointer to model region tags.
       - "stiffness variation": expression
     */
    pAttributeRefNode src_AttRefNode, snk_AttRefNode;
    amsi::cutPaste<pANode>(AttNode_childByType((pANode)cs,"stiffness source"),src_AttRefNode);
    amsi::cutPaste<pANode>(AttNode_childByType((pANode)cs,"stiffness sink"), snk_AttRefNode);

    pANode src_nd, snk_nk;
    amsi::cutPaste<pAttributeRefNode>(AttributeRefNode_value(src_AttRefNode),src_nd);
    amsi::cutPaste<pAttributeRefNode>(AttributeRefNode_value(snk_AttRefNode),snk_nd);

    std::vector<pModelItem> mdl_src_itms;
    std::vector<double> mdl_src_ents;
    amsi::getAssociatedModelItems(cs,src_nd,std::back_inserter(mdl_src_itms));
    std::transform(mdl_src_itms.begin(),mdl_src_itms.end(),std::back_inserter(mdl_src_ents),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    std::cout << "Stiffness sources affects model entities : ";
    for (auto mdl_ent = mdl_src_itms.begin(); mdl_ent != mdl_src_itms.end(); ++mdl_ent)
    {
      assert(ModelItem_isGentity(*mdl_ent));
      std::cout << ModelItem_tag(*mdl_ent) << " ";
    }
    std::cout << std::endl;

    std::vector<pModelItem> mdl_snk_itms;
    std::vector<double> mdl_snk_ents;
    amsi::getAssociatedModelItems(cs,snk_nd,std::back_inserter(mdl_snk_itms)); 
    std::transform(mdl_snk_itms.begin(),mdl_snk_itms.end(),std::back_inserter(mdl_snk_ents),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    std::cout << "Stiffness sources affects model entities : ";
    for (auto mdl_ent = mdl_snk_itms.begin(); mdl_ent != mdl_snk_itms.end(); ++mdl_ent)
    {
      assert(ModelItem_isGentity(*mdl_ent));
      std::cout << ModelItem_tag(*mdl_ent) << " ";
    }
    std::cout << std::endl;

    // Need to extract expression from "stiffness variation"
    stfvar = new StiffnessVariation(mdl_src_ents.begin(), mdl_src_ents.end(), mdl_snk_ents.begin(), mdl_snk_ents.end(),stf_vrtn_fld);
    return stfvar;
  }
  void StiffnessVariation::setStf_Vrtn_Fld()
  {
    apf::Mesh * msh = apf::getMesh(stf_vrtn_fld);
    int dm = msh->getDimension();
    for (auto snk = mdl_snk_ents.begin(); snk != mdl_snk_ents.end(); ++snk)
    {
      /* Identify faces that are adjacent to model region entity snk. */
      adj_mdl_snk_fcs GR_faces((pGRegion)snk);
      void *iter = 0;
      while (adj_mdl_snk_fc = PList_next(adj_mdl_snk_fcs, &iter))
      {
	for (auto src = mdl_src_ents.begin(); src != mdl_src_ents.end(); ++src)
	{
	  /* if adjacent model sink face is in closure of model source region. */
	  if (GR_inClosure((pGRegion)src,adj_mdl_snk_fc))
	  {
	    for (auto msh_rgn = amsi::beginClassified(msh,(pGRegion)src,dm); msh_rgn != amsi::endClassified(ent); ++msh_rgn)
	    {
	      (*this)->set_mdl_src_fc(adj_mdl_snk_fc);
	      apf::MeshElement * msh_elmt = apf::createMeshElement(msh,msh_rgn);
	      (*this)->process(msh_elmt);
	    }
	  }
	}
      }
    }
  }
  void StiffnessVariation::atPoint(apf::Vector3 const &p, double w, double dV)
  {
    /* Find global coordinates for Gauss Integration Point */
    apf::Mesh * msh = apf::getMesh(stf_vrtn_fld);
    int dm = msh->getDimension();
    apf::Vector3 gauss_pt;
    apf::getGaussPoint(apf::getType(msh),0,p,gauss_pt);
    /* Find closest point to model face entity that has been identified as a source. */
    double gp[3] = {0};
    double mdl_src_fc_xyz[3] = {0};
    double mdl_src_fc_uv[2] = {0};
    gauss_pt->toArray(gp);
    GF_closestPoint(mdl_src_fc,gp[],mdl_src_fc_xyz[],mdl_src_fc_uv[]);
    for (auto i=0; i<dm; ++i)
      dist += (gp[i] - mdl_fc_xyz[i]) * (gp[i] - mdl_fc_xyz[i]);
    dist = std::sqrt(dist);
    apf::setScalar(stf_vrtn,msh_ent,ip_integration_pt,f(dist));
    ip_itegration_pt++;
    //analysis->storeStiffnessVariation(me,dist);
  }
  
}
