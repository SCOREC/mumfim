#ifndef MUMFIM_STIFFNESS_VARIATION_H_
#define MUMFIM_STIFFNESS_VARIATION_H_
#include <apf.h>
#include <model_traits/CategoryNode.h>
#include <model_traits/ModelTrait.h>
#include <apfMesh.h>
#include <cstring>
#include <memory>
namespace mumfim
{
  class StiffnessVariation;
  std::unique_ptr<StiffnessVariation> buildStiffnessVariation(
      const mt::CategoryNode & nd,
      apf::Field * stf_vrtn_fld);
  class StiffnessVariation : public apf::Integrator
  {
    protected:
    apf::Field * stf_vrtn_fld;
    std::vector<apf::ModelEntity *> mdl_src_ents;
    std::vector<apf::ModelEntity *> mdl_snk_ents;
    apf::Vector3 target_xyz;
    apf::MeshEntity * msh_ent;
    int ip_integration_pt;
    std::shared_ptr<const mt::ScalarFunctionMT<1>> fn;
    public:
    StiffnessVariation(std::vector<apf::ModelEntity*> src_ents,
                       std::vector<apf::ModelEntity*> sink_ents,
                       apf::Field * stf_vrtn_fld,
                       std::shared_ptr<const mt::ScalarFunctionMT<1>> fn);
    void inElement(apf::MeshElement * msh_lmnt) final
    {
      msh_ent = apf::getMeshEntity(msh_lmnt);
      ip_integration_pt = 0;
    }
    void outElement() final {};
    void atPoint(apf::Vector3 const &, double, double) final;
    void set_target_xyz(const apf::Vector3 & xyz) { target_xyz = xyz; }
    void populate_stf_vrtn_fld();
    void calculate_mdl_fc_coord(apf::ModelEntity * mdl_fc, apf::Vector3& mdl_fc_xyz);
  };
}  // namespace mumfim
#endif
