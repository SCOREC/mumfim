namespace bio
{
  template <typename I>
    StiffnessVariation::StiffnessVariation(I mdl_src_ent_bgn, I mdl_src_ent_end, I mdl_snk_ent_bgn, I mdl_snk_end_end, apf::Field * f)
    , stf_vrtn_fld(n)
    , mdl_src_ents()
    , mdl_snk_ends()
    , src_mdl_fc(NULL)
    , ip_integration_pt(0)
  {
    std::copy(mdl_src_ent_bgn, mdl_src_ent_end, std::back_inserter(mdl_src_ents));
    std::copy(mdl_snk_ent_bgn, mdl_snk_ent_end, std::back_inserter(mdl_snk_ends));
  }
}
