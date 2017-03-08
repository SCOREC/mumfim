#ifndef BIO_VOLUME_CONVERGENCE_H_
#define BIO_VOLUME_CONVERGENCE_H_
#include "bioNonlinearTissue.h"
#include <simNonlinearAnalysis.h>
#include <apfMeasure.h>
namespace bio
{
  amsi::Convergence * buildBioConvergenceOperator(pACase ss, pANode cn, amsi::Iteration * it, apf::Field * fld);
  // Volume convergence class that considers accumulated volume of all regions where DeltaV = V-Vprev
  template <typename I>
  struct volCalc : public amsi::to_R1
  {
  public:
  private:
    double v0;
    double v;
    double pv;
    std::vector<apf::ModelEntity*> mdl_ents;
    apf::Field * u;
  };
  template <typename EPS>
    class VolumeConvergence : public amsi::UpdatingConvergence<EPS>
  {
  protected:
    int rnk;
    double v0;
    double v;
    double pv;
    amsi::Log vols;

    apf::Field * u;
    double calcVolume()
    {
      return amsi::measureDisplacedModelEntities(mdl_ents.begin(),mdl_ents.end(),u);
    }
  public:
    template <typename I>
      VolumeConvergence(I b, I n,EPS e,apf::Field * u_)
      : amsi::UpdatingConvergence<EPS>(e)
      , rnk(-1)
      , v(0.0)
      , pv(0.0)
      , vols(amsi::activateLog("volume"))
      , mdl_ents()
      , u(u_)
    {
      MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
      std::copy(b,n,std::back_inserter(mdl_ents));
      v = v0 = pv = calcVolume();
    }
    virtual void update()
    {
      pv = v;
      v = calcVolume();
      amsi::UpdatingConvergence<EPS>::update();
    }
    bool converged()
    {
      update();
      bool converged = false;
      double dv = abs(v - pv);
      // convergence based on volume change
      converged = dv < amsi::UpdatingConvergence<EPS>::eps * pv ;
      std::cout << "current volume: " << v << std::endl;
      std::cout << "previous volume: " << pv << std::endl;
      std::cout << "accumulated incremental volume convergence: " << std::endl
                << "\t" << dv << " < " << amsi::UpdatingConvergence<EPS>::eps * pv << std::endl
                << "\t" << (converged ? "TRUE" : "FALSE") << std::endl;
      return converged;
    }
    bool failed()
    {
      return false;
    }
    void log(int ldstp, int rnk)
    {
      if (rnk==0)
        amsi::log(vols) << ldstp << ", "
                        << "region tags..." << ", "
                        << v0 << ", "
                        << v << ", "
                        << v - pv << ", "
                        << v - v0 << std::endl;
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
