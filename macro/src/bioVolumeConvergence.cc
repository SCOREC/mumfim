#include "bioVolumeConvergence.h"
#include <amsiCasters.h>
#include <gmi.h>
namespace bio
{
  struct CalcDV : public amsi::to_R1
  {
    VolCalc * v;
    double operator()()
    {
      return fabs(v->getV() - v->getPV());
    }
  };
  struct CalcDV0 : public amsi::to_R1
  {
    VolCalc * v;
    double operator()()
    {
      return fabs(v->getV() - v->getV0());
    }
  };
  struct CalcDVPS : public amsi::to_R1
  {
    VolCalc * v;
    double operator()()
    {
      return fabs(v->getV() - v->getVPS());
    }
  };
  struct CalcPV : public amsi::to_R1
  {
    VolCalc * v;
    double operator()()
    {
      return v->getPV();
    }
  };
  struct CalcV : public amsi::to_R1
  {
    VolCalc * v;
    double operator()()
    {
      return v->getV();
    }
  };
  struct CalcV0 : public amsi::to_R1
  {
    VolCalc * v;
    double operator()()
    {
      return v->getV0();
    }
  };
  struct CalcVPS : public amsi::to_R1
  {
    VolCalc * v;
    double operator()()
    {
      return v->getVPS();
    }
  };
  amsi::Convergence * buildVolConvergenceOperator(pACase ss, pAttribute cn, amsi::Iteration * it, VolCalc * vl, apf::Field * fld)
  {
    pAttribute rgn_nd_att = Attribute_childByType(cn,"regions");
    pANode rgn_nd = AttributeRefNode_value((pAttributeRefNode)rgn_nd_att);
    std::vector<apf::ModelEntity*> mdl_ents;
    amsi::getAssociatedModelItems(ss,rgn_nd,std::back_inserter(mdl_ents));
    VolCalc cv(mdl_ents.begin(),mdl_ents.end(),fld);
    assert(cv == *vl);
    std::cout << "Volume convergence operator discovered, effects model entitites : ";
    gmi_model * mdl = apf::getMesh(fld)->getModel();
    for(auto mdl_ent = mdl_ents.begin(); mdl_ent != mdl_ents.end(); ++mdl_ent)
      std::cout << mdl->ops->tag(mdl,(gmi_ent*)*mdl_ent) << " ";
    std::cout << std::endl;
    pAttribute ref_att = Attribute_childByType(cn,"reference value");
    int ref_tp = AttributeInt_value((pAttributeInt)ref_att);
    amsi::to_R1 * dv = NULL;
    amsi::to_R1 * ref_v = NULL;
    if(ref_tp == 0) //  initial
    {
      CalcDV0 * v = new CalcDV0;
      v->v = vl;
      dv = v;
      CalcV0 * r = new CalcV0;
      r->v = vl;
      ref_v = r;
    }
    else if (ref_tp == 1) // load_step
    {
      CalcDVPS * v = new CalcDVPS;
      v->v = vl;
      dv = v;
      CalcVPS * r = new CalcVPS;
      r->v = vl;
      ref_v = r;
    }
    else if (ref_tp == 2) // iteration
    {
      CalcDV * v = new CalcDV;
      v->v = vl;
      dv = v;
      CalcVPS * r = new CalcVPS;
      r->v = vl;
      ref_v = r;
    }
    pAttribute eps_att = Attribute_childByType(cn,"epsilon");
    pAttribute cap_att = Attribute_childByType(cn,"iteration cap");
    amsi::SimUpdatingEpsilon * eps = new amsi::SimUpdatingEpsilon((pAttributeDouble)eps_att);
    if(cap_att)
      eps->setCap(AttributeInt_value((pAttributeInt)cap_att));
    return new amsi::UpdatingConvergence<amsi::to_R1*,decltype(eps),amsi::to_R1*>(it,dv,eps,ref_v);
  }
}
