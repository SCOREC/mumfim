#ifndef BIO_STIFFNESS_VARIATION_H_
#define BIO_STIFFNESS_VARIATION_H_

namespace bio
{
  class StiffnessVariation;
  StiffnessVariation * buildStiffnessVariation(pACase cs, pANode nd, apf::Numbering * un)
  class StiffnessVariation : public apf::integrator
  {
  protected:
    std::vector<apf::ModelEntity*> mdl_src_ents;
    std::vector<apf::ModelEntity*> mdl_snk_ents;
    apf::Numbering * nm;
    apf::Element * e;
    apf::MeshElement * me;
    int nen;
    int nedofs;
    virtual void _inElement(apf::MeshEntity * msh_ent) = 0;
  public:
    template <typename I>
      StiffnessVariation(I mdl_src_ent_bgn, I mdl_src_ent_end, I mdl_snk_ent_bgn, I mdl_snk_end_end, apf::Numbering * nm);
    void inElement(apf::MeshEntity * msh_ent);
    void outElement();
  };
  class Axon_StiffnessVariation : public StiffnessVariation
  {
  protected:
    apf::ModelEntity * src_mdl_fc;
    double Pt2fc_dist;
    virtual void _inElement(apf::MeshEntity *){};
  public:
    template <typename I>
      Axon_StiffnessVariation(I mdl_src_ent_bgn, I mdl_src_ent_end, I mdl_snk_ent_bgm, I mdl_snk_end_end, apf::Numbering * nm);
    void atPoint(apf::Vector3 const &, double, double);
    void set_src_mdl_fc(apf::ModelEntity * fc){src_mdl_fc = fc};
    double get_Pt2fc_dist(){return Pt2fc_dist};
  };
}
#include "bioStiffnessVariation_impl.h"
#endif
