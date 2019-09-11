#ifndef BIO_VOLUME_CONVERGENCE_H_
#define BIO_VOLUME_CONVERGENCE_H_
#include "bioNonlinearTissue.h"
#include <simNonlinearAnalysis.h>
#include <apfMeasure.h>
namespace bio
{
  class VolCalc : public amsi::Iteration, public amsi::PerStep
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
    virtual void iterate() 
    {
      pv = v;
      v = amsi::measureDisplacedModelEntities(mdl_ents.begin(),mdl_ents.end(),u);
      Iteration::iterate();
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
                                                  amsi::MultiIteration * it,
                                                  VolCalc * vl,
                                                  apf::Field * fld);
}
#endif
