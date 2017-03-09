#include "bioVolumeConvergence.h"
#include <amsiCasters.h>
#include <gmi.h>
namespace bio
{
  amsi::Convergence * buildBioConvergenceOperator(pACase ss, pANode cn, amsi::Iteration * it, apf::Field * fld)
  {
    pANode reg_nd_nd = AttNode_childByType(cn,"regions");
    pANode rgn_nd = AttInfoRefNode_value((pAttInfoRefNode)reg_nd_nd);
    std::vector<apf::ModelEntity*> mdl_ents;
    amsi::getAssociatedModelItems(ss,rgn_nd,std::back_inserter(mdl_ents));
    std::cout << "Volume convergence operator discovered, effects model entitites : ";
    gmi_model * mdl = apf::getMesh(fld)->getModel();
    for(auto mdl_ent = mdl_ents.begin(); mdl_ent != mdl_ents.end(); ++mdl_ent)
      std::cout << mdl->ops->tag(mdl,(gmi_ent*)*mdl_ent) << " ";
    std::cout << std::endl;
    pANode ref_nd = AttNode_childByType(cn,"reference");
    pANode eps_nd = AttNode_childByType(cn,"epsilon");
    int ref_tp = AttInfoInt_value((pAttInfoInt)ref_nd);
    VolCalc * dv = NULL;
    VolCalc * ref_v = NULL;
    if(ref_tp == 0) //  initial
    {
      dv = new CalcDV0(mdl_ents.begin(),mdl_ents.end(),fld);
      ref_v = new CalcV0(mdl_ents.begin(),mdl_ents.end(),fld);
    }
    else if (ref_tp == 1) // load_step
    {
      std::cerr << "ERROR: Cannot use load step as reference for volume convergence! (unimplemented)" << std::endl;
      //dv = new CalcDV0(mdl_ents.begin(),mdl_ents.end(),fld);
      //ref_v = new Calc0(mdl_ents.begin(),mdl_ents.end(),fld);
    }
    else if (ref_tp == 2) // iteration
    {
      dv = new CalcDV0(mdl_ents.begin(),mdl_ents.end(),fld);
      ref_v = new CalcPV(mdl_ents.begin(),mdl_ents.end(),fld);
    }
    amsi::SimUpdatingEpsilon * eps = new amsi::SimUpdatingEpsilon((pAttInfoDouble)eps_nd);
    return new amsi::UpdatingConvergence<VolCalc*,amsi::SimUpdatingEpsilon*,VolCalc*>(it,dv,eps,ref_v);
  }
}
