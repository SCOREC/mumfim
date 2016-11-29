#ifndef BIO_NONLINEARTISSUE_H_
#define BIO_NONLINEARTISSUE_H_
#include "bioLinearTissue.h"
#include <apfFEA.h>
#include "MicroFOMultiscaleTypes.h"
#include "RepresentVolElem.h" // should be able to take this out... (needed for RVE_Info struct)
#include "bioPreserveVolConstraintIntegrator.h"
#include "bioVolumeConstraintIntegrator.h"
#include <Solvers.h>
#include <NonLinElasticity.h>
#include <amsiMultiscale.h>
#include <amsiAnalysis.h>
#include <amsiUtil.h>
#include <apfsimWrapper.h>
#include <vector>
#include <list>
#include <string>
#include <iostream>
namespace bio
{
  class NonlinearTissue : public amsi::apfSimFEA
  {
  protected:
    std::vector<VolumeConstraintSurface*> vol_cnst;
    amsi::ElementalSystem * constitutive;
    apf::Field * delta_u;
  private:
    amsi::ElementalSystem * stress_strain_system;
    apf::Field * strs;
    apf::Field * rcvrd_strs;
    apf::Field * strn;
    apf::Field * previous_rve;
    double poisson_ratio;
    double youngs_modulus;
    std::vector<double> rgn_vols;
    std::vector<double> init_rgn_vols;
    std::vector<double> prev_rgn_vols;
    double dv_prev;
    int load_step;
    int iteration;
    int iteration_beta;
  public:
    NonlinearTissue(pGModel imdl,
                 pParMesh imsh,
                 pACase pd,
                 MPI_Comm cm = AMSI_COMM_SCALE);
    ~NonlinearTissue();
    virtual void ApplyBC_Dirichlet();
    virtual void RenumberDOFs();
    void getLoadOn(pGEntity ent, double * frc);
    void step();
    void iter() {iteration++;}
    apf::Numbering * getNumbering() { return apf_primary_numbering; }
    apf::Field * getdUField() { return delta_u; }
    apf::Field * getUField() {return apf_primary_field;}
    virtual void Assemble(amsi::LAS * las);
    virtual void UpdateDOFs(const double * );
    void UpdateLambda();
    void UpdateBeta(double);
    void setBeta(double);
    void setLambda(double);
    double getRgnVol(int idx){return rgn_vols[idx];}
    double getRgnVolInit(int idx){return init_rgn_vols[idx];}
    double getRgnVolPrev(int idx){return prev_rgn_vols[idx];}
    int getIteration(){return iteration;}
    void computeInitGuess(amsi::LAS * las)
    {
      LinearTissue lt(model,mesh,prob_def,analysis_comm);
      lt.setSimulationTime(T);
      LinearSolver(&lt,las);
      las->iter();
      apf::copyData(delta_u,lt.getField());
      apf::copyData(apf_primary_field,lt.getField());
      //amsi::PrintField(delta_u,std::cout).run();
    }
    void ComputeDispL2Norm(double &);
    void recoverSecondaryVariables(int);
    void storeStress(apf::MeshElement * me, double * stress);
    void storeStrain(apf::MeshElement * me, double * strain);
    double getPoissonRatio(){return poisson_ratio;}
    double getYoungsModulus(){return youngs_modulus;}
    int numVolumeConstraints()
    {
      return vol_cnst.size();
    }
    VolumeConstraintSurface * getVolumeConstraint(int idx)
    {
      if((int)vol_cnst.size() >= idx)
        return vol_cnst[idx];
      else
        return NULL;
    }
    void updateVolumes()
    {
      int ii = 0;
      double vol_glb = 0.0;
      double init_vol_glb = 0.0;
      for(auto it = vol_cnst.begin(); it != vol_cnst.end(); ++it)
      {
        rgn_vols[ii] = amsi::measureVol_pGFace((*it)->getFace(),(*it)->getRegionTag(),part,apf_primary_field);
//      rgn_vols[ii] = amsi::measureDisplacedEntity((*it)->getRegion(),part,apf_primary_field);
        vol_glb += rgn_vols[ii];
        init_vol_glb += init_rgn_vols[ii];
        ii++;
      }
      std::cout<<"current volume (updateVolumes()) = "<<vol_glb<<std::endl;
      std::cout<<"initial volume (updateVolumes()) = "<<init_vol_glb<<std::endl;
      for(auto it = vol_cnst.begin(); it != vol_cnst.end(); ++it)
        (*it)->setVol(vol_glb);
    }
    void updatePrevVolumes()
    {
      prev_rgn_vols = rgn_vols;
      double prev_vol_glb = 0.0;
      for (uint ii = 0; ii < prev_rgn_vols.size(); ii++)
        prev_vol_glb += prev_rgn_vols[ii];
      std::cout<<"previous volume (updatePrevVolumes()) = "<<prev_vol_glb<<std::endl;
      for (auto it = vol_cnst.begin(); it != vol_cnst.end(); ++it)
        (*it)->setPrevVol(prev_vol_glb);
    }
    void resetLambda()
    {
      for(std::vector<VolumeConstraintSurface*>::iterator cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
        (*cnst)->setLambda(0.0);
    }
    void updateConstraints();
    void updateConstraintsAccm();
    void updateConstraintsAccm_Incrmt();
    /// record parameters for imposing Volume constraint.
    void logCnstrntParams(int ldstp, int iteration, int rnk);
  protected:
    void calcInitVolumes()
    {
      int ii = 0;
      init_rgn_vols.resize(vol_cnst.size());
      for(auto it = vol_cnst.begin(); it != vol_cnst.end(); ++it)
      {
        init_rgn_vols[ii] = amsi::measureVol_pGFace((*it)->getFace(),(*it)->getRegionTag(),part,apf_primary_field);
        ii++;
      }
    }
  };
  // Volume convergence class that considers accumulated volume of all regions where DeltaV = V-Vprev
  class VolumeConvergenceAccm_Incrmt : public amsi::Convergence
  {
  protected:
    NonlinearTissue * ts;
    double eps;
    amsi::Log vols;
  public:
    VolumeConvergenceAccm_Incrmt(NonlinearTissue * tssu, double e)
      : amsi::Convergence()
      , ts(tssu)
      , eps(e)
      , vols(amsi::activateLog("volume"))
    { }
    bool converged()
    {
      ts->updateVolumes();
      int rgns = ts->numVolumeConstraints();
      bool converged = true;
      double vi = 0.0; // Accumulated Current Volume
      double vp = 0.0; // Accumulated Previous Volume
      double dv = 0.0; // Accumulated Volume Difference
      if (rgns == 0)
        vp = 1.0; // so that we do not run into divide by zero error.
      for(int ii = 0; ii < rgns; ii++)
      {
        // prev_rgn_vol initialized to init_rgn_vol,
        // therefore, dv_rgn = vi_rgn - initial v_rgn at ldstp0, iteration 0.
        double vp_rgn = ts->getRgnVolPrev(ii);
        double vi_rgn = ts->getRgnVol(ii);
        double dv_rgn = vi_rgn - vp_rgn;
        vi += vi_rgn;
        vp += vp_rgn;
        dv += dv_rgn;
      }
      // convergence based on volume change
      converged = std::abs(dv) < eps * vp  ;
      std::cout << "current volume: " << vi << std::endl;
      std::cout << "previous volume: " << vp << std::endl;
      std::cout << "accumulated incremental volume convergence: " << std::endl
                << "\t" << dv << " < " << eps * vp << std::endl
                << "\t" << (converged ? "TRUE" : "FALSE") << std::endl;
      if(!converged)
        ts->updateConstraintsAccm_Incrmt();
      return converged;
    }
    bool failed()
    {
      return false;
    }
    void log(int ldstp, int iteration, int rnk)
    {
      int rgns = ts->numVolumeConstraints();
      double v0 = 0;
      double vi = 0;
      double vp = 0; ///<Accumulated Previous Volume.
      for (int ii = 0; ii < rgns; ii++)
      {
        double vi_rgn = ts->getRgnVol(ii);
        double v0_rgn = ts->getRgnVolInit(ii);
        double vp_rgn = ts->getRgnVolPrev(ii);
        v0 += v0_rgn;
        vi += vi_rgn;
        vp += vp_rgn;
      }
      if (rnk==0)
        amsi::log(vols) << ldstp << ", "
                        << iteration << ", "
                        << "entire domain" << ", "
                        << v0 << ", "
                        << vi << ", "
                        << vi - vp << ", "
                        << vi - v0 << std::endl;
    }
  };
}
#endif
