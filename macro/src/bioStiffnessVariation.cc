#include "bioStiffnessVariation.h"
#include <amsiCasters.h>
#include <gmi.h>
#include <simClassified.h>
#include <math.h> // std::sqrt()
namespace bio
{
  StiffnessVariation * buildStiffnessVariation(pACase cs, pANode nd, apf::Field * stf_vrtn_fld)
  {
    StiffnessVariation * stfvar = NULL;

    /* cs = "problem definition" in simmodeler.
       nd = "stiffness gradient" in simmodeler. nd has three children.
       - "source"   : of type pAttributeRefNode
       - "sink"     : of type pAttributeRefNode
       - "function" : expression
     */
    pAttributeRefNode src_AttRefNode, snk_AttRefNode;
    src_AttRefNode = (pAttributeRefNode)AttNode_childByType(nd,"source");
    snk_AttRefNode = (pAttributeRefNode)AttNode_childByType(nd,"sink");

    pANode src_nd, snk_nd;
    src_nd = AttributeRefNode_value(src_AttRefNode);
    snk_nd = AttributeRefNode_value(snk_AttRefNode);

    std::vector<pModelItem> mdl_src_itms;
    std::vector<apf::ModelEntity*> mdl_src_ents;
    amsi::getAssociatedModelItems(cs,src_nd,std::back_inserter(mdl_src_itms));
    std::transform(mdl_src_itms.begin(),mdl_src_itms.end(),std::back_inserter(mdl_src_ents),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    /* New implementation
       std::vector<apf::ModelEntity*> mdl_src_ents;
       amsi::getAssociatedModelItems(cs,src_nd,std::back_inserter(mdl_src_ents));
     */
    std::cout << "Stiffness sources include model entities : ";
    for (auto mdl_ent = mdl_src_itms.begin(); mdl_ent != mdl_src_itms.end(); ++mdl_ent)
    {
      assert(ModelItem_isGEntity(*mdl_ent));
      std::cout << ModelItem_tag(*mdl_ent) << " ";
    }
    std::cout << std::endl;

    std::vector<pModelItem> mdl_snk_itms;
    std::vector<apf::ModelEntity*> mdl_snk_ents;
    amsi::getAssociatedModelItems(cs,snk_nd,std::back_inserter(mdl_snk_itms)); 
    std::transform(mdl_snk_itms.begin(),mdl_snk_itms.end(),std::back_inserter(mdl_snk_ents),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    /* New implementation
       std::vector<apf::ModelEntity*> mdl_snk_ents;
       amsi::getAssociatedModelItems(cs,snk_nd,std::back_inserter(mdl_snk_ents));
     */
    std::cout << "Stiffness sources include model entities : ";
    for (auto mdl_ent = mdl_src_itms.begin(); mdl_ent != mdl_src_itms.end(); ++mdl_ent)
    {
      assert(ModelItem_isGEntity(*mdl_ent));
      std::cout << ModelItem_tag(*mdl_ent) << " ";
    }
    std::cout << std::endl;

    pANode fn_nd;
    fn_nd = AttNode_childByType(nd,"function");
    stfvar = new StiffnessVariation(mdl_src_ents.begin(), mdl_src_ents.end(), mdl_snk_ents.begin(), mdl_snk_ents.end(), stf_vrtn_fld, fn_nd);
    return stfvar;
  }
  void StiffnessVariation::populate_stf_vrtn_fld()
  {
    /* Outline
       1. Iterate through model sink regions (std::vector<apf::ModelEntity*> mdl_snk_ents).
       2. For each region of model sink regions, identify adjacent faces (gmi_set * adj_mdl_snk_fcs).
       3. For each adjacent face (gmi_ent * mdl_snk_fc), iterate through model source regions (std::vector<apf::ModelEntity*> mdl_src_ents).
       4. If adjacent face is in closure of source region in model source regions, iterate through mesh regions that are classified on sink region.
       5. For each classified mesh region, set the target model face to adjacent face and process MeshElement.
     */
    apf::Mesh * msh = apf::getMesh(stf_vrtn_fld);
    int dm = msh->getDimension();
    gmi_model * mdl = msh->getModel();
    for (auto snk = mdl_snk_ents.begin(); snk != mdl_snk_ents.end(); ++snk)
    {
      /* Identify faces that are adjacent to model region entity snk. snk if of type apf::ModelEntity**. */
      gmi_set * adj_mdl_snk_fcs = mdl->ops->adjacent(mdl,(gmi_ent*)(*snk),2);
      for (int ii = 0; ii < adj_mdl_snk_fcs->n; ++ii)
      {
	gmi_ent * mdl_snk_fc = adj_mdl_snk_fcs->e[ii];
	for (auto mdl_src = mdl_src_ents.begin(); mdl_src != mdl_src_ents.end(); ++mdl_src)
	{
	  /* if adjacent model sink face is in closure of model source region.
	     - mdl        : type gmi_model*.
	     - mdl_src    : type apf::ModelEntity **.
	     - mdl_snk_fc : type gmi_ent *.
	   */
	  if (gmi_is_in_closure_of(mdl,(gmi_ent*)(*mdl_src),mdl_snk_fc))
	  {
	    for (auto msh_rgn = amsi::beginClassified(msh,*mdl_src,dm); msh_rgn != amsi::endClassified(msh_rgn); ++msh_rgn)
	    {
	      this->set_mdl_src_fc((apf::ModelEntity*)mdl_snk_fc);
	      apf::MeshElement * msh_elmt = apf::createMeshElement(msh,*msh_rgn);
	      this->process(msh_elmt);
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
    gmi_model * mdl = msh->getModel();
    apf::Vector3 gauss_pt;
    apf::getGaussPoint(msh->getType(msh_elmt),0,ip_integration_pt,gauss_pt);
    /* Find closest point to model face entity that has been identified as a source. */
    double gp[3] = {0};
    double mdl_src_fc_xyz[3] = {0};
    double mdl_src_fc_uv[2] = {0};
    gauss_pt.toArray(gp);
    gmi_closest_point(mdl,(gmi_ent*)mdl_src_fc,gp,mdl_src_fc_xyz,mdl_src_fc_uv);
    double dist = 0.0;
    for (auto i=0; i<dm; ++i)
      dist += (gp[i] - mdl_src_fc_xyz[i]) * (gp[i] - mdl_src_fc_xyz[i]);
    dist = std::sqrt(dist);
    // The expression from pANode fn is assumed to have the form 1.0 - dist/x, where x is specified in simmodeler.
    double x_val = 10.0; //assume variation distance is 10 microns.    
    apf::setScalar(stf_vrtn_fld,msh_elmt,ip_integration_pt,1.0-dist/x_val);
    ip_integration_pt++;
  }
  
}
