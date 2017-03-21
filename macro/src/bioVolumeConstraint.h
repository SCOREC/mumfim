#ifndef BIO_VOLUME_CONSTRAINT_H_
#define BIO_VOLUME_CONSTRAINT_H_
#include "bioConstraint.h"
#include <amsiDofHolder.h>
#include <amsiLAS.h>
#include <ElementalSystem.h>
#include <simAnalysis.h>
#include <apfShape.h>
#include <apfSIM.h>
#include <cstring>
namespace bio
{
  class VolumeConstraint;
  VolumeConstraint * buildVolumeConstraint(pACase cs, pANode nd, apf::Numbering * un);
  class VolumeConstraint : public Constraint, public apf::Integrator
  {
  protected:
    std::vector<apf::ModelEntity*> mdl_ents;
    apf::Numbering * nm;
    apf::Element * e;
    apf::MeshElement * me;
    int nen;
    int nedofs;
    double vol;
    double init_vol;
    double prev_vol;
    virtual void _update() = 0;
    virtual void _inElement(apf::MeshEntity * me) = 0;
  public:
    template <typename I>
      VolumeConstraint(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm);
    virtual double calcVolume();
    virtual void apply(amsi::LAS * las) = 0;
    void update();
    void inElement(apf::MeshEntity * me);
    void outElement();
    double getVolume() { return vol; }
    double getInitVolume() { return init_vol; }
    double getPrevVolume() { return prev_vol; }
    bool includesBodyForces() { return true; }
  };
  class LagrangeConstraint_Volume : public VolumeConstraint , public amsi::DofHolder
  {
  protected:
    apf::DynamicMatrix dVdu;
    apf::DynamicMatrix d2Vdu2;
    double lambda; // Lagrange multiplier for volume Preservation.
    double beta;   // Penalty parameter for Augmented Lagrangian Method.
    virtual void _inElement(apf::MeshEntity * me);
    virtual void _update();
  public:
    template <typename I>
      LagrangeConstraint_Volume(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double b);
    void apply(amsi::LAS * las);
    void atPoint(apf::Vector3 const &p, double w, double dV);
    apf::DynamicMatrix& getdVdu(){return dVdu;}
    apf::DynamicMatrix& getd2Vdu2(){return d2Vdu2;}
    double getLambda() { return lambda; }
  };
  class PenaltyConstraint_VolumeSurface : public VolumeConstraint
  {
  protected:
    apf::Mesh * msh;
    apf::ModelEntity * crt_rgn;
    apf::ModelEntity * crt_fc;
    apf::DynamicMatrix dVdu;
    double lambda;
    double beta;
    virtual void _update() {}
    virtual void _inElement(apf::MeshEntity *);
    void calcdVdu(apf::Vector3 const & pt0, apf::Vector3 const & pt1, apf::Vector3 const & pt2);
    //void calcd2Vdu2(apf::Vector3 const & pt0, apf::Vector3 const & pt1, apf::Vector3 const & pt2);
  public:
    template <typename I>
      PenaltyConstraint_VolumeSurface(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double b);
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
      PenaltyConstraint_Volume(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double b);
    void apply(amsi::LAS *) {}
    void atPoint(apf::Vector3 const &, double, double) {}
  };
  class LagrangeConstraint_VolumeSurface : public VolumeConstraint
  {
  protected:
    double beta;
    virtual void _update();
  public:
    template <typename I>
      LagrangeConstraint_VolumeSurface(I mdl_ent_bgn, I mdl_ent_end, apf::Numbering * nm, double b);
    void apply(amsi::LAS *, apf::Numbering *) {}
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
    LagrangeVolumeConstraint_Surface(apf::ModelEntity * me, apf::Numbering * nm);
    void apply(amsi::LAS * las, apf::Mesh * msh, apf::Numbering * nm);
    void atPoint(apf::Vector3 const &p, double w, double dV);
    void calcG(apf::Vector3 const & pt0,
               apf::Vector3 const & pt1,
               apf::Vector3 const & pt2,
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
}
#include "bioVolumeConstraint_impl.h"
#endif
