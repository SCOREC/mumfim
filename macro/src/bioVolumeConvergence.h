#ifndef BIO_VOLUME_CONVERGENCE_H_
#define BIO_VOLUME_CONVERGENCE_H_
#include "bioNonlinearTissue.h"
#include <simNonlinearAnalysis.h>
#include <apfMeasure.h>
namespace bio
{
  class MultiscaleConvergence : public amsi::MultiConvergence
  {
  private:
    amsi::ControlService * cs;
    size_t cplg;
  public:
    template <typename I>
      MultiscaleConvergence(I bgn, I end, size_t c)
      : amsi::MultiConvergence(bgn,end)
      , cs(amsi::ControlService::Instance())
      , cplg(c)
    { }
    virtual bool converged()
    {
      bool rslt = MultiscaleConvergence::converged();
      cs->scaleBroadcast(cplg,&rslt);
      return rslt;
    }
  };
  class VolCalc : public amsi::PerIter, public amsi::PerStep
  {
  public:
    template <typename I>
      VolCalc(I bgn, I end, apf::Field * _u)
      : v0(0.0)
      , vps(0.0)
      , v(0.0)
      , pv(0.0)
      , mdl_ents()
      , u(_u)
    {
      std::copy(bgn,end,std::back_inserter(mdl_ents)); // use amsi caster
      v = vps = pv = v0 = amsi::measureDisplacedModelEntities(mdl_ents.begin(),mdl_ents.end(),u);
    }
    void iter()
    {
      pv = v;
      v = amsi::measureDisplacedModelEntities(mdl_ents.begin(),mdl_ents.end(),u);
    }
    void step()
    {
      vps = amsi::measureDisplacedModelEntities(mdl_ents.begin(),mdl_ents.end(),u);
    }
    double getV0() const
    {
      return v0;
    }
    double getVPS() const
    {
      return vps;
    }
    double getV() const
    {
      return v;
    }
    double getPV() const
    {
      return pv;
    }
    bool operator==(VolCalc & o)
    {
      return std::equal(mdl_ents.begin(),mdl_ents.end(),o.mdl_ents.begin()) && u == o.u;
    }
    protected:
    double v0;
    double vps;
    double v;
    double pv;
    std::vector<apf::ModelEntity*> mdl_ents;
    apf::Field * u;
  };
  amsi::Convergence * buildVolConvergenceOperator(pACase ss,
                                                  pAttribute cn,
                                                  amsi::Iteration * it,
                                                  VolCalc * vl,
                                                  apf::Field * fld);
}
#endif
