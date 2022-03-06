#include "bioStiffnessVariation.h"
#include <amsiCasters.h>
#include <apfWrapper.h>
#include <gmi.h>
#include <math.h>  // std::sqrt()
#include <model_traits/CategoryNode.h>
#include <pumi.h>
#include <memory>
#include "bioModelTraits.h"
namespace mumfim
{
  // pass in the "stiffness gradient category node
  std::unique_ptr<StiffnessVariation> buildStiffnessVariation(
      const mt::CategoryNode & nd,
      apf::Field * stf_vrtn_fld)
  {
    const auto * src_cat = nd.FindCategoryNodeByType("sources");
    const auto * sink_cat = nd.FindCategoryNodeByType("sinks");
    if (src_cat == nullptr || sink_cat == nullptr)
    {
      std::cerr << "stiffness variation requires a \"sources\" and \"sinks\" "
                   "category.\n";
      exit(1);
    }
    src_cat = src_cat->FindCategoryNodeByType("stiffness source");
    sink_cat->FindCategoryNodeByType("stiffness sink");
    if (src_cat == nullptr || sink_cat == nullptr)
    {
      std::cerr
          << "sources category must contain a \"stiffness source\"  and \"stiffness sink\" category.\n";
      exit(1);
    }

    std::vector<apf::ModelEntity*> mdl_src_ents;
    for(const auto& src: src_cat->GetModelTraitNodes())
    {
      GetModelTraitNodeGeometry(apf::getMesh(stf_vrtn_fld), &src, mdl_src_ents);
    }
    std::cout << "Stiffness sources include model entities : ";
    for (auto* mdl_ent : mdl_src_ents)
    {
      std::cout << apf::getMesh(stf_vrtn_fld)->getModelTag(mdl_ent) << " ";
    }
    std::cout << std::endl;

    std::vector<apf::ModelEntity*> mdl_snk_ents;
    for(const auto& sink: sink_cat->GetModelTraitNodes())
    {
      GetModelTraitNodeGeometry(apf::getMesh(stf_vrtn_fld), &sink,
                                mdl_src_ents);
    }
    std::cout << "Stiffness sinks include model entities : ";
    for (auto* mdl_ent : mdl_snk_ents)
    {

      std::cout << apf::getMesh(stf_vrtn_fld)->getModelTag(mdl_ent) << " ";
    }
    std::cout << std::endl;
    const auto * fn_nd = nd.FindModelTraitNode("function");
    if (fn_nd == nullptr)
    {
      std::cerr << "stiffness gradient must have a \"function model trait\".\n";
      exit(1);
    }
    auto fn = std::dynamic_pointer_cast<const mt::ScalarFunctionMT<1>>(
        fn_nd->GetModelTraits()[0].second);
    if (fn == nullptr)
    {
      std::cerr
          << "stiffness gradient function model trait must be a function that "
             "takes a single variable (the distance from the source)\n";
      exit(1);
    }
    return std::make_unique<StiffnessVariation>(
        std::move(mdl_src_ents), std::move(mdl_snk_ents), stf_vrtn_fld, fn);
  }
  void StiffnessVariation::populate_stf_vrtn_fld()
  {
    /* Outline
       1. Iterate through model regions labeled as "src"
       (std::vector<apf::ModelEntity*> mdl_src_ents).
       2. Identify model faces that are adjacent to each model region in
       mdl_src_ents (gmi_set * adj_mdl_src_fcs).
       3. For each adjacent face (gmi_ent * mdl_src_fc), calculate average
       coordinate.
       4. Iterate through model regions labeled as "snk"
       (std::vector<apf::ModelEntity*> mdl_snk_ents).
       5. If mdl_src_fc is in closure of snk.
       6. Iterate through msh regions classified on snk and process each msh
       region.
    */
    /* 1. Iterate through model regions labeled as "src" */
    apf::Mesh * msh = apf::getMesh(stf_vrtn_fld);
    int dm = msh->getDimension();
    gmi_model * mdl = msh->getModel();
    apf::MeshEntity * msh_ent;
    auto * it = msh->begin(dm);
    while ((msh_ent = msh->iterate(it)))
    {
      for (auto * src : mdl_src_ents)
      {
        /* 2. identify model faces that are adjacent to each model region in
         * mdl_src_ents. */
        gmi_set * adj_mdl_src_fcs =
            mdl->ops->adjacent(mdl, (gmi_ent *)(src), 2);
        // for each face on the source region
        for (int ii = 0; ii < adj_mdl_src_fcs->n; ++ii)
        {
          gmi_ent * mdl_src_fc =
              adj_mdl_src_fcs
                  ->e[ii];  // extract ii^th adj. face from adj_mdl_src_fcs.
          apf::Vector3 fc_xyz;
          /* 3. Calculate average coordinate of adjacent face, mdl_src_fc */
          calculate_mdl_fc_coord((apf::ModelEntity *)mdl_src_fc, fc_xyz);
          /* 4. Iterate through model regions labeled as "snk" */
          // for (auto snk = mdl_snk_ents.begin(); snk != mdl_snk_ents.end();
          // ++snk)
          for (auto * snk : mdl_snk_ents)
          {
            /* 5. if adjacent model src face is in closure of model sink region,
               process the model sink region
               - mdl        : type gmi_model*.
               - mdl_src    : type apf::ModelEntity **.
               - mdl_snk_fc : type gmi_ent *.
            */
            std::cout << "face " << gmi_tag(mdl, mdl_src_fc)
                      << " of src mdl rgn " << gmi_tag(mdl, (gmi_ent *)(src))
                      << " in closure of snk mdl rgn "
                      << gmi_tag(mdl, (gmi_ent *)(snk)) << "?" << std::endl;
            // if mdl_src_face in closure of the sink region
            if (gmi_is_in_closure_of(mdl, mdl_src_fc, (gmi_ent *)(snk)))
            {
              std::cout << " Yes, coord. of face " << gmi_tag(mdl, mdl_src_fc)
                        << " is " << fc_xyz[0] << "," << fc_xyz[1] << ","
                        << fc_xyz[2] << std::endl;
              /* 6. Iterate through msh regions that are classified in "snk" and
               * process each msh region. */
              if (msh->toModel(msh_ent) == snk)
              {
                set_target_xyz(fc_xyz);
                apf::MeshElement * msh_elmt =
                    apf::createMeshElement(msh, msh_ent);
                inElement(msh_elmt);
                process(msh_elmt);
                apf::destroyMeshElement(msh_elmt);
              }
            }
          }
        }
      }
    }
    msh->end(it);
  }
  void StiffnessVariation::calculate_mdl_fc_coord(apf::ModelEntity * mdl_fc,
                                                  apf::Vector3 & mdl_fc_xyz)
  {
    apf::Mesh * mesh = apf::getMesh(stf_vrtn_fld);
    std::set<gmi_ent *> face_verts;
    // get the centroid of the model face
    mesh->snapToModel(mdl_fc, {0.5, 0.5, 0}, mdl_fc_xyz);
  }
  void StiffnessVariation::atPoint(apf::Vector3 const & p, double, double)
  {
    /* Find global coordinates for Gauss Integration Point */
    apf::Mesh * msh = apf::getMesh(stf_vrtn_fld);
    int dm = msh->getDimension();
    apf::Vector3 gauss_pt_xyz;
    apf::mapLocalToGlobal(apf::createMeshElement(msh, msh_ent), p,
                          gauss_pt_xyz);
    double gp[3] = {0};
    gauss_pt_xyz.toArray(gp);
    double dist = 0.0;
    for (auto i = 0; i < dm; ++i)
      dist += (gp[i] - target_xyz[i]) * (gp[i] - target_xyz[i]);
    dist = std::sqrt(dist);
    // The expression from pANode fn is assumed to have the form 1-dist/x, where
    // x is specified in simmodeler.
    double val = (*fn)(dist);
    // src region is processed multiple times, therefore take maximum value of
    // 1-dist/x_val among all process instances.
    double stf_vrtn_coeff =
        std::max(val, apf::getScalar(stf_vrtn_fld, msh_ent, ip_integration_pt));
    apf::setScalar(stf_vrtn_fld, msh_ent, ip_integration_pt, stf_vrtn_coeff);
    ip_integration_pt++;
  }
  StiffnessVariation::StiffnessVariation(std::vector<apf::ModelEntity*> src_ents,
                     std::vector<apf::ModelEntity*> sink_ents,
                     apf::Field * stf_vrtn_fld,
                     std::shared_ptr<const mt::ScalarFunctionMT<1>> fn)
      : apf::Integrator(apf::getShape(stf_vrtn_fld)->getOrder())
      , stf_vrtn_fld(stf_vrtn_fld)
      , mdl_src_ents(std::move(src_ents))
      , mdl_snk_ents(std::move(sink_ents))
      , target_xyz()
      , ip_integration_pt(0)
      , fn(std::move(fn))
  {
  }
}  // namespace mumfim
