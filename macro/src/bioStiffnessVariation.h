#ifndef BIO_STIFFNESS_VARIATION_H_
#define BIO_STIFFNESS_VARIATION_H_
#include <simAnalysis.h>
#include <apfSIM.h>
#include <cstring>
#include <apf.h>
namespace bio
{
  class StiffnessVariation;
  StiffnessVariation * buildStiffnessVariation(pACase cs, pANode nd, apf::Field * stf_vrtn);
  class StiffnessVariation : public apf::Integrator
  {
  protected:
    apf::Field * stf_vrtn_fld;
    std::vector<apf::ModelEntity*> mdl_src_ents;
    std::vector<apf::ModelEntity*> mdl_snk_ents;
    double * target_xyz;
    apf::MeshEntity * msh_ent;
    int ip_integration_pt;
    pANode fn;
  public:
    template <typename I>
      StiffnessVariation(I mdl_src_ent_bgn, I mdl_src_ent_end, I mdl_snk_ent_bgn, I mdl_snk_end_end, apf::Field * stf_vrtn_fld, pANode fn);
    void inElement(apf::MeshEntity * msh_ent_in){msh_ent = msh_ent_in; ip_integration_pt=0;}
    void outElement(){};
    void atPoint(apf::Vector3 const &, double, double);
    void set_target_xyz(double * xyz){target_xyz = xyz;}
    void populate_stf_vrtn_fld();
    void calculate_mdl_fc_coord(apf::ModelEntity * mdl_fc, double * mdl_fc_xyz);
  };
}
#include "bioStiffnessVariation_impl.h"
#endif
