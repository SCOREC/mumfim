#ifndef BIO_VOLUME_CONSTRAINT_H_
#define BIO_VOLUME_CONSTRAINT_H_
#include <amsiLAS.h>
#include <ElementalSystem.h>
#include <simAnalysis.h>
#include <apfShape.h>
#include <apfSIM.h>
#include <cstring>
namespace bio
{
  class VolumeConstraint;
  VolumeConstraint * buildVolumeConstraint(pACase cs, pANode nd, apf::Field * u);
  class VolumeConstraint : public amsi::ElementalSystem
  {
  protected:
    std::vector<apf::ModelEntity*> mdl_ents;
    double vol;
    double init_vol;
    double prev_vol;
    virtual void _inElement(apf::MeshEntity * me);
  public:
    template <typename I>
      VolumeConstraint(I mdl_ent_bgn, I mdl_ent_end, apf::Field * fld);
    virtual double calcVolume();
    virtual void update();
    virtual void _update() = 0;
    virtual void apply(amsi::LAS * las, apf::Numbering * nm) = 0;
    double getVolume() { return vol; }
    double getInitVolume() { return init_vol; }
    double getPrevVolume() { return prev_vol; }
    bool includesBodyForces() { return true; }
  };
  class LagrangeConstraint_Volume : public VolumeConstraint
  {
  protected:
    apf::DynamicMatrix dVdu;
    apf::DynamicMatrix d2Vdu2;
    double dV;
    double lambda; // Lagrange multiplier for volume Preservation.
    double beta;   // Penalty parameter for Augmented Lagrangian Method.
    void _inElement(apf::MeshEntity * me);
  public:
    template <typename I>
      LagrangeConstraint_Volume(I mdl_ent_bgn, I mdl_ent_end, apf::Field * fld);
    virtual void update(); // called each iteration..
    virtual void apply(amsi::LAS * las, apf::Mesh * msh, apf::Numbering * nm);
    void atPoint(apf::Vector3 const &p, double w, double dV);
    apf::DynamicMatrix& getdVdu(){return dVdu;}
    apf::DynamicMatrix& getd2Vdu2(){return d2Vdu2;}
    double getLambda() { return lambda; }
    double getdV() { return dV; }
  };
  class PenaltyConstraint_Volume : public VolumeConstraint
  {
  public:
    template <typename I>
      PenaltyConstraint_Volume(I mdl_ent_bgn, I mdl_ent_end, apf::Field * fld, int o);
  };
  class LagrangeConstraint_VolumeSurface : public VolumeConstraint
  {

  };
  class PenaltyConstraint_VolumeSurface : public VolumeConstraint
  {

  };
  /*
  class VolumeConstraintADMM : public VolumeConstraint
  {
  protected:
    std::vector<double> prev_vols;
    int elem_num; // keeps count of element number within region.
  public:
    VolumeConstraintADMM(apf::ModelEntity * me, apf::Field * fld, int o);
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
    LagrangeVolumeConstraint_Surface(apf::ModelEntity * me, apf::Field * fld, int o);
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
  class PenaltyVolumeConstraint_Surface : public LagrangeVolumeConstraint_Surface
  {
  public:
    PenaltyVolumeConstraint_Surface(apf::ModelEntity * me, apf::Field * fld, int o)
      : LagrangeVolumeConstraint_Surface(me,fld,o)
    { }
  };
  */
}
#endif
