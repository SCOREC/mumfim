#ifndef MUMFIM_TISSUEMULTISCALEANALYSIS_H_
#define MUMFIM_TISSUEMULTISCALEANALYSIS_H_
#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiPETScLAS.h>
#include <apf.h>
#include <model_traits/CategoryNode.h>
#include <stdexcept>
#include <string>
#include "LinearTissue.h"
#include "MultiscaleTissue.h"
#include "TissueAnalysis.h"
#include "VolumeConvergence.h"
namespace mumfim
{
  class MultiscaleTissueIteration : public amsi::Iteration
  {
    protected:
    MultiscaleTissue * tssu;
    amsi::LAS * las;
    Iteration * fem_iter;

    public:
    MultiscaleTissueIteration(MultiscaleTissue * a, amsi::LAS * l)
        : tssu(a), las(l), fem_iter(amsi::buildLinearFEMIteration(a, l))
    {
    }
    ~MultiscaleTissueIteration() { delete fem_iter; }
    void iterate();
  };
  class MultiscaleTissueAnalysis : public TissueAnalysis
  {
    public:
    MultiscaleTissueAnalysis(apf::Mesh * mesh,
                             std::unique_ptr<mt::CategoryNode> analysis_case,
                             MPI_Comm cm);
    virtual void init();
    virtual void run();
    virtual void finalizeStep();

    private:
    size_t cplng;
  };
}  // namespace mumfim
#endif
