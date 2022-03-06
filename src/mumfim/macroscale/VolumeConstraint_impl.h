#include <apfNumbering.h>
namespace mumfim
{
  template <typename I>
    VolumeConstraint::VolumeConstraint(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * n)
    : apf::Integrator(apf::getShape(apf::getField(n))->getOrder())
    , lg(amsi::activateLog("constraints"))
    , mdl_ents()
    , nm(n)
    , nen(0)
    , nedofs(0)
    , vol(0.0)
    , prev_vol(0.0)
    , load_vol(0.0)
    , init_vol(0.0)
  {
    std::copy(mdl_ent_bgn,mdl_ent_end,std::back_inserter(mdl_ents));
    vol = init_vol = prev_vol = calcVolume();
  }
  template <typename I>
    LagrangeConstraint_Volume::LagrangeConstraint_Volume(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double b)
    : VolumeConstraint(mdl_ent_bgn, mdl_ent_end, nm)
    , dVdu()
    , d2Vdu2()
    , lambda(0.0)
    , beta(b)
  { }
  template <typename I>
    VolumeSurfaceConstraint::VolumeSurfaceConstraint(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm)
    : VolumeConstraint(mdl_ent_bgn,mdl_ent_end,nm)
    , msh(apf::getMesh(apf::getField(nm)))
    , crt_rgn()
    , crt_fc()
  { }
  template <typename I>
    LagrangeConstraint_VolumeSurface::LagrangeConstraint_VolumeSurface(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double l, double b)
    : VolumeSurfaceConstraint(mdl_ent_bgn,mdl_ent_end,nm)
    , lambda(l)
    , beta(b)
  { }
  template <typename I>
    PenaltyConstraint_VolumeSurface::PenaltyConstraint_VolumeSurface(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double b)
    : VolumeSurfaceConstraint(mdl_ent_bgn, mdl_ent_end, nm)
    , beta(b)
  { }
  /*
  template <typename I>
    PenaltyConstraint_Volume(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double b)
    : VolumeConstraint(mdl_ent_bgn, mdl_ent_end, nm)
    , beta(b)
  { }
  template <typename I>
    LagrangeConstraint_VolumeSurface(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double b)
    : VolumeConstraint(mdl_ent_bgn, mdl_ent_end, nm)
    , beta(b)
  { }
  */
}
