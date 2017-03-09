#include "bioVolumeConvergence.h"
#include <amsiCasters.h>
namespace bio
{
  amsi::Convergence * buildBioConvergenceOperator(pACase ss, pANode cn, Iteration * it, apf::Field * fld)
  {
    amsi::Convergence * cnvg = NULL;
    pANode rgn_nd = AttInfoRefNode_value(cn);
    std::vector<apf::ModelEntity*> mdl_ents;
    std::vector<pModelItem> mdl_itms;
    amsi::getAssociatedModelItems(ss,rgn_nd,std::back_inserter(mdl_itms));
    std::transform(mdl_itms.begin(),mdl_itms.end(),std::back_inserter(mdl_ents),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
    std::cout << "Volume convergence operator discovered, effects model entitites : ";
    for(auto mdl_ent = mdl_itms.begin(); mdl_ent != mdl_itms.end(); ++mdl_ent)
    {
      assert(ModelItem_isGEntity(*mdl_ent));
      std::cout << ModelItem_tag(*mdl_ent) << " ";
    }
    std::cout << std::endl;
    pANode type = AttNode_childByType(cn,"type");
    pANode ref  = AttNode_childByType(cn,"reference");
    pANode eps = AttNode_childByType(cn,"epsilon");
    assert(type);
    assert(ref);
    assert(eps);
    //todo: parse types
    VolCalc * dv = new CalcDV(mdl_itms.begin(),mdl_itms.end(),u);
    amsi::SimUpdatingEpsiloneps = new amsi::SimUpdatingEpsilon((pAttInfoDouble)eps,it);
    VolCalc * vp = new CalcPC(mdl_itms.begin(),mdl_itms.end(),u);
    return new VolumeConvergence<VolCalc*,amsi::SimUpdatingEpsilon*,VolCalc*>(dv,eps,vp);
  }
}
