#ifndef BIO_STIFFNESS_VARIATION_H_
#define BIO_STIFFNESS_VARIATION_H_

namespace bio
{
  class StiffnessVariation;
  StiffnessVariation * buildStiffnessVariation(pACase cs, pANode nd, apf::Field * stf_vrtn)
  class StiffnessVariation : public apf::integrator
  {
  protected:
    apf::Field * stf_vrtn_fld;
    std::vector<apf::ModelEntity*> mdl_src_ents;
    std::vector<apf::ModelEntity*> mdl_snk_ents;
    apf::ModelEntity * src_mdl_fc;
    int ip_integration_pt;
  public:
    template <typename I>
      StiffnessVariation(I mdl_src_ent_bgn, I mdl_src_ent_end, I mdl_snk_ent_bgn, I mdl_snk_end_end, apf::Field * stf_vrtn);
    void inElement(apf::MeshEntity * msh_ent){ip_integration_pt=0};
    void outElement(){};
    void atPoint(apf::Vector3 const &, double, double);
    void set_mdl_src_fc(apf::ModelEntity * fc){src_mdl_fc = fc};
  };
}
#include "bioStiffnessVariation_impl.h"
#endif
