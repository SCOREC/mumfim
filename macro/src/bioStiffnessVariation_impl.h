namespace bio
{
  template <typename I>
    StiffnessVariation::StiffnessVariation(I mdl_src_ent_bgn, I mdl_src_ent_end, I mdl_snk_ent_bgn, I mdl_snk_ent_end, apf::Field * f, pANode fn)
    : apf::Integrator(apf::getShape(f)->getOrder())
    , stf_vrtn_fld(f)
    , mdl_src_ents()
    , mdl_snk_ents()
    , mdl_fc_xyz()
    , ip_integration_pt(0)
    , fn(fn)
  {
    std::copy(mdl_src_ent_bgn, mdl_src_ent_end, std::back_inserter(mdl_src_ents));
    std::copy(mdl_snk_ent_bgn, mdl_snk_ent_end, std::back_inserter(mdl_snk_ents));
  }
}
