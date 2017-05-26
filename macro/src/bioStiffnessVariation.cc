#include "bioStiffnessVariation.h"
#include <amsiCasters.h>
#include <gmi.h>
#include <simClassified.h>
#include <simWrapper.h>
#include <apfWrapper.h>
#include <math.h> // std::sqrt()
namespace bio
{
  StiffnessVariation * buildStiffnessVariation(pACase cs, pANode nd, apf::Field * stf_vrtn_fld)
  {
    StiffnessVariation * stfvar = NULL;
    /* cs = "problem definition" in simmodeler.
       nd = "stiffness gradient" in simmodeler. nd has three children.
       - "source"   : of type pAttInfoRefNode
       - "sink"     : of type pAttInfoRefNode
       - "function" : expression
    */
    amsi::describeNode(nd);
    pANode src_AttRefNode, snk_AttRefNode;
    src_AttRefNode = AttNode_childByType(nd,"sources");
    snk_AttRefNode = AttNode_childByType(nd,"sinks");
    amsi::describeNode(src_AttRefNode);
    pANode src_nd, snk_nd;
    src_nd = AttInfoRefNode_value((pAttInfoRefNode)src_AttRefNode);
    snk_nd = AttInfoRefNode_value((pAttInfoRefNode)snk_AttRefNode);
    amsi::describeNode(src_nd);
    std::vector<apf::ModelEntity*> mdl_src_ents;
    amsi::getAssociatedModelItems(cs,src_nd,std::back_inserter(mdl_src_ents));
    std::cout << "Stiffness sources include model entities : ";
    gmi_model * mdl = apf::getMesh(stf_vrtn_fld)->getModel();
    for (auto mdl_ent = mdl_src_ents.begin(); mdl_ent != mdl_src_ents.end(); ++mdl_ent)
      std::cout << mdl->ops->tag(mdl,(gmi_ent*)*mdl_ent) << " ";
    std::cout << std::endl;
    std::vector<apf::ModelEntity*> mdl_snk_ents;
    amsi::getAssociatedModelItems(cs,snk_nd,std::back_inserter(mdl_snk_ents));
    std::cout << "Stiffness sinks include model entities : ";
    for (auto mdl_ent = mdl_snk_ents.begin(); mdl_ent != mdl_snk_ents.end(); ++mdl_ent)
      std::cout << mdl->ops->tag(mdl,(gmi_ent*)*mdl_ent) << " ";
    std::cout << std::endl;
    pAttribute stf_att = AttCase_attrib(cs,"stiffness variation");
    pAttributeDouble fn = (pAttributeDouble)Attribute_childByType(stf_att,"function");
    //fn_nd = AttNode_childByType(nd,"function");
    stfvar = new StiffnessVariation(mdl_src_ents.begin(), mdl_src_ents.end(), mdl_snk_ents.begin(), mdl_snk_ents.end(), stf_vrtn_fld, fn);
    return stfvar;
  }
  void StiffnessVariation::populate_stf_vrtn_fld()
  {
    /* Outline
       1. Iterate through model regions labeled as "src" (std::vector<apf::ModelEntity*> mdl_src_ents).
       2. Identify model faces that are adjacent to each model region in mdl_src_ents (gmi_set * adj_mdl_src_fcs).
       3. For each adjacent face (gmi_ent * mdl_src_fc), calculate average coordinate.
       4. Iterate through model regions labeled as "snk" (std::vector<apf::ModelEntity*> mdl_snk_ents).
       5. If mdl_src_fc is in closure of snk.
       6. Iterate through msh regions classified on snk and process each msh region.
    */
    /* 1. Iterate through model regions labeled as "src" */
    apf::Mesh * msh = apf::getMesh(stf_vrtn_fld);
    int dm = msh->getDimension();
    gmi_model * mdl = msh->getModel();
    for (auto src = mdl_src_ents.begin(); src != mdl_src_ents.end(); ++src)
    {
      /* 2. identify model faces that are adjacent to each model region in mdl_src_ents. */
      gmi_set * adj_mdl_src_fcs = mdl->ops->adjacent(mdl,(gmi_ent*)(*src),2);
      for (int ii = 0; ii < adj_mdl_src_fcs->n; ++ii)
      {
        gmi_ent * mdl_src_fc = adj_mdl_src_fcs->e[ii]; // extract ii^th adj. face from adj_mdl_src_fcs.
        double fc_xyz[3] = {0.0};
        /* 3. Calculate average coordinate of adjacent face, mdl_src_fc */
        calculate_mdl_fc_coord((apf::ModelEntity*)mdl_src_fc, fc_xyz);
        /* 4. Iterate through model regions labeled as "snk" */
        for (auto snk = mdl_snk_ents.begin(); snk != mdl_snk_ents.end(); ++snk)
        {
          /* 5. if adjacent model src face is in closure of model sink region, process the model sink region
             - mdl        : type gmi_model*.
             - mdl_src    : type apf::ModelEntity **.
             - mdl_snk_fc : type gmi_ent *.
          */
          std::cout<<"face "<<gmi_tag(mdl,mdl_src_fc)<<" of src mdl rgn "<<gmi_tag(mdl,(gmi_ent*)(*src))<< " in closure of snk mdl rgn "
                   <<gmi_tag(mdl,(gmi_ent*)(*snk))<<"?"<<std::endl;
          if (gmi_is_in_closure_of(mdl,mdl_src_fc,(gmi_ent*)(*snk)))
          {
            std::cout<<" Yes, coord. of face "<<gmi_tag(mdl,mdl_src_fc)<<" is "<<fc_xyz[0]<<","<<fc_xyz[1]<<","<<fc_xyz[2]<<std::endl;
            /* 6. Iterate through msh regions that are classified in "snk" and process each msh region. */
            for (auto msh_rgn = amsi::beginClassified(msh,*snk,dm); msh_rgn != amsi::endClassified(msh_rgn); ++msh_rgn)
            {
              this->set_target_xyz(fc_xyz);
              this->inElement(*msh_rgn);
              apf::MeshElement * msh_elmt = apf::createMeshElement(msh,*msh_rgn);
              this->process(msh_elmt);
            }
          }
        }
      }
    }
  }
  void StiffnessVariation::calculate_mdl_fc_coord(apf::ModelEntity * mdl_fc, double * mdl_fc_xyz)
  {
    apf::Mesh * msh = apf::getMesh(stf_vrtn_fld);
    auto bgn = amsi::beginClassified(msh, mdl_fc, 0);
    auto end = amsi::endClassified(bgn);
    amsi::getAvgFieldValue(msh->getCoordinateField(),bgn,end,mdl_fc_xyz);
  }
  void StiffnessVariation::atPoint(apf::Vector3 const &p, double w, double dV)
  {
    /* Find global coordinates for Gauss Integration Point */
    apf::Mesh * msh = apf::getMesh(stf_vrtn_fld);
    int dm = msh->getDimension();
    apf::Vector3 gauss_pt_xyz;
    apf::mapLocalToGlobal(apf::createMeshElement(msh,msh_ent),p,gauss_pt_xyz);
    double gp[3] = {0};
    gauss_pt_xyz.toArray(gp);
    double dist = 0.0;
    for (auto i=0; i<dm; ++i)
      dist += (gp[i] - target_xyz[i]) * (gp[i] - target_xyz[i]);
    dist = std::sqrt(dist);
    // The expression from pANode fn is assumed to have the form 1-dist/x, where x is specified in simmodeler.
    double val = AttributeDouble_evalDT(fn,dist);
    // src region is processed multiple times, therefore take maximum value of 1-dist/x_val among all process instances.
    double  stf_vrtn_coeff = std::max(val,apf::getScalar(stf_vrtn_fld,msh_ent,ip_integration_pt));
    apf::setScalar(stf_vrtn_fld,msh_ent,ip_integration_pt,stf_vrtn_coeff);
    ip_integration_pt++;
  }
}
