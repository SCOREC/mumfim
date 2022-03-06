#ifndef BIO_VOLUME_CONSTRAINT_H_
#define BIO_VOLUME_CONSTRAINT_H_
#include <ElementalSystem.h>
#include <amsiDofHolder.h>
#include <amsiLAS.h>
#include <amsiNonlinearAnalysis.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <model_traits/AssociatedModelTraits.h>
#include <cstring>
#include <memory>
#include "bioConstraint.h"
#include <amsiReporter.h>
namespace bio
{
  class VolumeConstraint;
  std::unique_ptr<VolumeConstraint> buildVolumeConstraint(
      const mt::CategoryNode & mt,
      apf::Numbering * un);
  class VolumeConstraint : public apf::Integrator,
                           public amsi::PerIter,
                           public amsi::PerStep
  {
    protected:
    amsi::Log lg;
    std::vector<apf::ModelEntity *> mdl_ents;
    apf::Numbering * nm;
    apf::Element * e;
    apf::MeshElement * me;
    int nen;
    int nedofs;
    double vol;
    double prev_vol;
    double load_vol;
    double init_vol;
    virtual void _iter() = 0;
    virtual void _step() = 0;
    virtual void _inElement(apf::MeshElement * me) = 0;

    public:
    template <typename I>
    VolumeConstraint(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm);
    virtual double calcVolume();
    virtual void apply(amsi::LAS * las) = 0;
    void iter();
    void step();
    void inElement(apf::MeshElement * me);
    void outElement();
    double getVolume() { return vol; }
    double getInitVolume() { return init_vol; }
    double getPrevVolume() { return prev_vol; }
    bool includesBodyForces() { return true; }
    virtual ~VolumeConstraint() = default;
  };
  class LagrangeConstraint_Volume : public VolumeConstraint,
                                    public amsi::DofHolder
  {
    protected:
    apf::DynamicMatrix dVdu;
    apf::DynamicMatrix d2Vdu2;
    double lambda;  // Lagrange multiplier for volume Preservation.
    double beta;    // Penalty parameter for Augmented Lagrangian Method.
    virtual void _inElement(apf::MeshElement * me);
    virtual void _iter();
    virtual void _step();

    public:
    template <typename I>
    LagrangeConstraint_Volume(I mdl_ent_bgn,
                              I mdl_ent_end,
                              apf::Numbering * nm,
                              double b);
    void apply(amsi::LAS * las);
    void atPoint(apf::Vector3 const & p, double w, double dV);
    apf::DynamicMatrix & getdVdu() { return dVdu; }
    apf::DynamicMatrix & getd2Vdu2() { return d2Vdu2; }
  };
  // class PenaltyConstraint_Volume : public VolumeConstraint , public
  // amsi::DofHolder
  class VolumeSurfaceConstraint : public VolumeConstraint
  {
    protected:
    apf::Mesh * msh;
    apf::ModelEntity * crt_rgn;
    apf::ModelEntity * crt_fc;

    public:
    template <typename I>
    VolumeSurfaceConstraint(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm);
  };
  class LagrangeConstraint_VolumeSurface : public VolumeSurfaceConstraint
  {
    protected:
    apf::DynamicMatrix dVdu;
    double lambda;
    double beta;
    virtual void _iter();
    virtual void _step(){};
    virtual void _inElement(apf::MeshElement *);

    public:
    template <typename I>
    LagrangeConstraint_VolumeSurface(I mdl_ent_bgn,
                                     I mdl_ent_end,
                                     apf::Numbering * nm,
                                     double l,
                                     double b);
    void apply(amsi::LAS *);
    void atPoint(apf::Vector3 const &, double, double);
  };
  class PenaltyConstraint_VolumeSurface : public VolumeSurfaceConstraint
  {
    protected:
    apf::DynamicMatrix dVdu;
    double beta;
    virtual void _iter() {}
    virtual void _step() {}
    virtual void _inElement(apf::MeshElement *);

    public:
    template <typename I>
    PenaltyConstraint_VolumeSurface(I mdl_ent_bgn,
                                    I mdl_ent_end,
                                    apf::Numbering * nm,
                                    double b);
    void apply(amsi::LAS *);
    void atPoint(apf::Vector3 const &, double, double);
  };
  /*
  class PenaltyConstraint_Volume : public VolumeConstraint
  {
  protected:
    double beta;
    virtual void _update() {};
  public:
    template <typename I>
      PenaltyConstraint_Volume(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering *
  nm, double b); void apply(amsi::LAS *) {} void atPoint(apf::Vector3 const &,
  double, double) {}
  };
  class LagrangeConstraint_VolumeSurface : public VolumeConstraint
  {
  protected:
    double beta;
    virtual void _update();
  public:
    template <typename I>
      LagrangeConstraint_VolumeSurface(I mdl_ent_bgn, I mdl_ent_end,
  apf::Numbering * nm, double b); void apply(amsi::LAS *, apf::Numbering *) {}
    void atPoint(apf::Vector3 const &, double, double) {}
  };
  */
  /*
  class VolumeConstraintADMM : public VolumeConstraint
  {
  protected:
    std::vector<double> prev_vols;
    int elem_num; // keeps count of element number within region.
  public:
    VolumeConstraintADMM(apf::ModelEntity * me, apf::Numbering ** fld, int o);
    void apply(amsi::LAS * las, apf::Mesh * msh, apf::Numbering * nm);
    void atPoint(apf::Vector3 const &p, double w, double dV);
    //void setPrevVol(int idx, double vol) { prev_vols[idx] = vol; }
    //void getPrevVols(std::vector<double> & v) { v = prev_vol; }
  };
  */
  /**
   * VolumeConstraintSurface class that contains implementation
   *  of volume constraint based on surface mesh. Calculations
   *  are based on Hong et al., "Fast Volume Preservation for a
   *  Mass-Spring System," IEEE Computer Graphics and Applications (2006)
   */
  /*
  class LagrangeVolumeConstraint_Surface : public VolumeConstraint
  {
  protected:
    bool should_update_G;
    bool should_update_H;
    double prev_vol;
    int elem_num; // keeps count of element number within region.
    virtual void _inElement(apf::MeshElement * me);
  public:
    LagrangeVolumeConstraint_Surface(apf::ModelEntity * me, apf::Numbering *
  nm); void apply(amsi::LAS * las, apf::Mesh * msh, apf::Numbering * nm); void
  atPoint(apf::Vector3 const &p, double w, double dV); void calcG(apf::Vector3
  const & pt0, apf::Vector3 const & pt1, apf::Vector3 const & pt2,
               apf::DynamicMatrix & dVdu);
    void calcH(apf::Vector3 const & pt0,
               apf::Vector3 const & pt1,
               apf::Vector3 const & pt2,
               apf::DynamicMatrix & dVdu);
    void update();
    void setLambda(double l) { lambda = l; }
    void updateG(bool uf) { should_update_G = uf; }
    void updateH(bool uf) { should_update_H = uf; }
    void setVol(double v) { vol = v; }
    void setInitVol(double vol) { init_vol = vol;}
    void setPrevVol(double vol) { prev_vol = vol;}
  };
  */
}  // namespace bio
#include "bioVolumeConstraint_impl.h"
#endif
