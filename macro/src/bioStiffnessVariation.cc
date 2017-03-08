#include "bioStiffnessVariation.h"
namespace bio
{
  StiffnessVariation * buildStiffnessVariation(pAcase cs, pANode nd, apf::Numbering * un)
  {
    StiffnessVariation * stfvar = NULL;
    
    /* The input pANode nd contains a attribute node of "stiffness variation" type. */
    std::vector<pANode> mdl_src_nds;
    std::vector<pANode> mdl_snk_nds;
    amsi::cutPaste<pANode>(AttNode_childByType(nd,"stiffness source"),std::back_inserter(mdl_src_nds));
    amsi::cutPaste<pANode>(AttNode_childByType(nd,"stiffness sink"),std::back_inserter(mdl_snk_nds));

    std::vector<pModelItem> mdl_src_itms;
    std::vector<pModelItem> mdl_snk_itms;
    std::vector<double> mdl_src_ents;
    std::vector<double> mdl_snk_ents;
    amsi::getAssociatedModelItems(cs,mdl_src_nds,std::back_inserter(mdl_src_itms));
    std::transform(mdl_src_itms.begin(),mdl_src_itms.end(),std::back_inserter(mdl_src_ents),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    std::cout << "Stiffness sources affects model entities : ";
    for (auto mdl_ent = mdl_src_itms.begin(); mdl_ent != mdl_src_itms.end(); ++mdl_ent)
    {
      assert(ModelItem_isGentity(*mdl_ent));
      std::cout << ModelItem_tag(*mdl_ent) << " ";
    }
    std::cout << std::endl;
    
    amsi::getAssociatedModelItems(cs,mdl_snk_nds,std::back_inserter(mdl_snk_itms)); 
    std::transform(mdl_snk_itms.begin(),mdl_snk_itms.end(),std::back_inserter(mdl_snk_ents),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    std::cout << "Stiffness sources affects model entities : ";
    for (auto mdl_ent = mdl_snk_itms.begin(); mdl_ent != mdl_snk_itms.end(); ++mdl_ent)
    {
      assert(ModelItem_isGentity(*mdl_ent));
      std::cout << ModelItem_tag(*mdl_ent) << " ";
    }
    std::cout << std::endl;
    stfvar = new Axon_StiffnessVariation(mdl_src_ents.begin(), mdl_src_ents.end(), mdl_snk_ents.begin(), mdl_snk_ents.end(),un);
    return stfvar;
  }
  void StiffnessVariation::inElement(apf::MeshEntity * msh_ent)
  {
    apf::Field * fld = apf::getField(nm);
    me = apf::createMeshElement(apf::getMesh(fld),msh_ent);
    e = apf::createElement(fld,me);
    nedofs = nen * apf::countComponents(fld);
    _inElement(m);
  }
  void StiffnessVariation::outElement()
  {
    apf::destroyMeshElement(me);
    apf::destroyElement(e);
  }
  double StiffnessVariation::calcClosestPt()
  {
    double dist = 0.0;
    for (auto 
  }
  void Axon_StiffnessVariation::atPoint(apf::Vector3 const &p, double w, double dV)
  {
    /* Find global coordinates for Gauss Integration Point */
    apf::Mesh * msh = apf::getMesh(apf::getField(nm));
    int dm = msh->getDimension();
    apf::Vector3 gauss_pt;
    apf::getGaussPoint(apf::getType(msh),0,p,gauss_pt);
    /* Find closest point to model face entity that has been identified as a source. */
    double gp[3] = {0};
    double mdl_fc_src_xyz[3] = {0};
    double mdl_fc_src_uv[2] = {0};
    gauss_pt->toArray(gp);
    GF_closestPoint(mdl_fc_src,gp[],mdl_fc_src_xyz[],mdl_fc_src_uv[]);
    for (auto i=0; i<dm; ++i)
      dist += (gp[i] - mdl_fc_xyz[i]) * (gp[i] - mdl_fc_xyz[i]);
    dist = std::sqrt(dist);
  }
  
}
