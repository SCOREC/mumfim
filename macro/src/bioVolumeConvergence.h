#ifndef BIO_VOLUME_CONVERGENCE_H_
#define BIO_VOLUME_CONVERGENCE_H_
#include "bioNonlinearTissue.h"
#include <simNonlinearAnalysis.h>
#include <apfMeasure.h>
namespace bio
{
  amsi::Convergence * buildBioConvergenceOperator(pACase ss, pANode cn, amsi::Iteration * it, apf::Field * fld);
  struct VolCalc : public amsi::to_R1
  {
  public:
    template <typename I>
      VolCalc(I bgn, I end, apf::Field * _u)
      : v0(0.0)
      , v(0.0)
      , pv(0.0)
      , mdl_ents()
      , u(_u)
    {
      std::copy(bgn,end,std::back_inserter(mdl_ents));
      v = pv = v0 = amsi::measureDisplacedModelEntities(mdl_ents.begin(),mdl_ents.end(),u);
    }
    void update()
    {
      pv = v;
      v = amsi::measureDisplacedModelEntities(mdl_ents.begin(),mdl_ents.end(),u);
    }
    protected:
    double v0;
    double v;
    double pv;
    std::vector<apf::ModelEntity*> mdl_ents;
    apf::Field * u;
  };
  struct CalcDV : public VolCalc
  {
    double operator()()
    {
      update();
      return abs(v - pv);
    }
  };
  struct CalcPV : public VolCalc
  {
    double operator()()
    {
      update();
      return pv;
    }
  };
  struct CalcV : public VolCalc
  {
    double operator()()
    {
      update();
      return v;
    }
  };
  struct CalcV0 : public VolCalc
  {
    double operator()()
    {
      return v0;
    }
  };
  /*
  /// Volume convergence class that considers current volume - initial volume of each individual region
  class VolumeConvergence : public amsi::Convergence
  {
  protected:
    MultiscaleTissue * ts;
    double eps;
    amsi::Log vols;
  public:
    VolumeConvergence(MultiscaleTissue * tssu, double e)
      : amsi::Convergence()
      , ts(tssu)
      , eps(e)
      , vols(amsi::activateLog("volume"))
    { }
    bool converged()
    {
      int rgns = ts->numVolumeConstraints();
      bool converged = true;
      for(int ii = 0; ii < rgns; ii++)
      {
//      VolumeConstraintADMM * cnstrnt = ts->getVolumeConstraint(ii);
        VolumeConstraintSurface * cnstrnt = ts->getVolumeConstraint(ii);
        double v0 = amsi::measureEntity(cnstrnt->getFace(),ts->getPart(),ts->getMesh());
        double vi = amsi::measureDisplacedEntity(cnstrnt->getFace(),ts->getPart(),ts->getUField());
        double dv = vi - v0;
        // convergence based on volume change
        converged = std::abs(dv) < eps  * v0;
        std::cout << "incremental volume " << ii << " convergence: " << std::endl
                  << "\t" << dv << " < " << eps * v0 << std::endl
                  << "\t" << (converged ? "TRUE" : "FALSE") << std::endl;
        if(!converged)
          break;
      }
      return converged;
    }
    bool failed()
    {
      return false;
    }
    void log(int ldstp, int iteration, int rnk)
    {
      int rgns = ts->numVolumeConstraints();
      for (int ii = 0; ii < rgns; ii++)
      {
//        VolumeConstraintADMM * cnstrnt = ts->getVolumeConstraint(ii);
        VolumeConstraintSurface * cnstrnt = ts->getVolumeConstraint(ii);
        double v0 = amsi::measureEntity(cnstrnt->getFace(),ts->getPart(),ts->getMesh());
        double vi = amsi::measureDisplacedEntity(cnstrnt->getFace(),ts->getPart(),ts->getUField());
        if (rnk==0)
          amsi::log(vols) << ldstp << ", "
                          << iteration << ", "
                          << GEN_tag(cnstrnt->getFace()) << ", "
                          << v0 << ", "
                          << vi << ", "
                          << vi - v0 << std::endl;
      }
    }
  };
  /// Volume convergence class that considers accumulated volume of all regions
  class VolumeConvergenceAccm : public amsi::Convergence
  {
  protected:
    MultiscaleTissue * ts;
    double eps;
    amsi::Log vols;
  public:
    VolumeConvergenceAccm(MultiscaleTissue * tssu, double e)
      : amsi::Convergence()
      , ts(tssu)
      , eps(e)
      , vols(amsi::activateLog("volume"))
    { }
    bool converged()
    {
      int rgns = ts->numVolumeConstraints();
      bool converged = true;
      double vi = 0.0; ///<Accumulated Current Volume
      double v0 = 0.0; ///<Accumulated Initial Volume
      double dv = 0.0; ///<Accumulated Volume Difference
      // If volume constraint is turned off, volume convergence criteria always true.
      if (rgns == 0)
        v0 = 1.0;
      for(int ii = 0; ii < rgns; ii++)
      {
        double v0_rgn = ts->getRgnVolInit(ii);
        double vi_rgn = ts->getRgnVol(ii);
        double dv_rgn = vi_rgn - v0_rgn;
        vi += vi_rgn;
        v0 += v0_rgn;
        dv += dv_rgn;
      }
      // convergence based on volume change
      converged = std::abs(dv) < eps * v0;
      std::cout << "current volume: " << vi << std::endl;
      std::cout << "initial volume: " << v0 << std::endl;
      std::cout << "accumulated volume convergence: " << std::endl
                << "\t" << dv << " < " << eps * v0  << std::endl
                << "\t" << (converged ? "TRUE" : "FALSE") << std::endl;
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
      for (int ii = 0; ii < rgns; ii++)
      {
        double v0_rgn = ts->getRgnVolInit(ii);
        double vi_rgn = ts->getRgnVol(ii);
        v0 += v0_rgn;
        vi += vi_rgn;
      }
      if (rnk==0)
        amsi::log(vols) << ldstp << ", "
                        << iteration << ", "
                        << "entire domain" << ", "
                        << v0 << ", "
                        << vi << ", "
                        << vi - v0 << std::endl;
    }
  };
  */
}
#endif
