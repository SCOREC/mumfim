#include "bioVolumeConvergence.h"
#include <amsiCasters.h>
#include <gmi.h>
namespace bio
{
  amsi::Convergence * buildBioConvergenceOperator(pACase ss, pAttribute cn, amsi::Iteration * it, apf::Field * fld)
  {
    pAttribute rgn_nd_att = Attribute_childByType(cn,"regions");
    pANode rgn_nd = AttributeRefNode_value((pAttributeRefNode)rgn_nd_att);
    std::vector<apf::ModelEntity*> mdl_ents;
    amsi::getAssociatedModelItems(ss,rgn_nd,std::back_inserter(mdl_ents));
    std::cout << "Volume convergence operator discovered, effects model entitites : ";
    gmi_model * mdl = apf::getMesh(fld)->getModel();
    for(auto mdl_ent = mdl_ents.begin(); mdl_ent != mdl_ents.end(); ++mdl_ent)
      std::cout << mdl->ops->tag(mdl,(gmi_ent*)*mdl_ent) << " ";
    std::cout << std::endl;
    pAttribute ref_att = Attribute_childByType(cn,"reference value");
    int ref_tp = AttributeInt_value((pAttributeInt)ref_att);
    VolCalc * dv = NULL;
    VolCalc * ref_v = NULL;
    if(ref_tp == 0) //  initial
    {
      dv = new CalcDV0(mdl_ents.begin(),mdl_ents.end(),fld);
      ref_v = new CalcV0(mdl_ents.begin(),mdl_ents.end(),fld);
    }
    else if (ref_tp == 1) // load_step
    {
      dv = new CalcDVPS(mdl_ents.begin(),mdl_ents.end(),fld);
      ref_v = new CalcVPS(mdl_ents.begin(),mdl_ents.end(),fld);
    }
    else if (ref_tp == 2) // iteration
    {
      dv = new CalcDV(mdl_ents.begin(),mdl_ents.end(),fld);
      ref_v = new CalcPV(mdl_ents.begin(),mdl_ents.end(),fld);
    }
    pAttribute eps_att = Attribute_childByType(cn,"epsilon");
    return new amsi::UpdatingConvergence<VolCalc*,amsi::SimUpdatingEpsilon*,VolCalc*>(it,
                                                                                      dv,
                                                                                      new amsi::SimUpdatingEpsilon((pAttributeDouble)eps_att),
                                                                                      ref_v);
  }
}
