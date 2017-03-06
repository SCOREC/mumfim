#include "bioVolumeConvergence.h"
namespace bio
{
  amsi::Convergence * buildBioConvergenceOperator(pACase ss, pANode cn, const int & it)
  {
    amsi::Convergence * cnvg = NULL;
    char * tp = AttNode_imageClass(cn);
    if(std::string("volume convergence").compare(tp) == 0)
    {
      std::vector<apf::ModelEntity*> mdl_ents;
      std::vector<pModelItem> mdl_itms;
      amsi::getAssociatedModelItems(ss,cn,std::back_inserter(mdl_itms));
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
      amsi::UpdatingEpsilon(eps,it)
        VolumeConvergence<amsi::UpdatingEpsilon>(
    }
    /*
    else
      cnfg = amsi::buildConvergenceOperator(ss,cn);
    */
    Sim_deleteString(tp);
  }
}
