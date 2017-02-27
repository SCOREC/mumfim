namespace bio
{
  template <typename I>
    VolumeConstraint::VolumeConstraint(I mdl_ent_bgn,
                                       I mdl_ent_end,
                                       apf::Field * fld)
    : ElementalSystem(fld,fld->getShape()->getOrder())
    , mdl_ents()
    , vol(0.0)
    , init_vol(0.0)
    , prev_vol(0.0)
  {
    std::copy(mdl_ent_bgn,mdl_ent_end,std::back_inserter(mdl_ents));
    vol = init_vol = prev_vol = calcVolume();
  }
  template <typename I>
    LagrangeConstraint_Volume(I mdl_ent_bgn,
                              I mdl_ent_end,
                              apf::Field * fld)
    : VolumeConstraint(mdl_ent_bgn,
                       mdl_ent_end,
                       fld)
    , dVdu()
    , d2Vdu2()
    , dV(0.0)
    , lambda(0.0)
    , beta(0.0)
  { }
}
