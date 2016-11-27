#ifndef BIO_TISSUE_ANALYSIS_H_
#define BIO_TISSUE_ANALYSIS_H_
#include "bioNonlinearTissue.h"
#include <amsiNonlinearAnalysis.h>
namespace bio
{
  class TissueIteration : public amsi::Iteration
  {
  protected:
    NonlinearTissue * tssu;
    amsi::LAS * las;
  public:
    TissueIteration(NonlinearTissue * t, amsi::LAS * l)
      : amsi::Iteration()
      , tssu(t)
      , las(l)
    {}
    virtual void iterate()
    {
      LinearSolver(tssu,las);
      las->iter();
    }
  };
}
#endif
