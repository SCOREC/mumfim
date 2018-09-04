#ifndef BIO_TISSUEMULTISCALEANALYSIS_H_
#define BIO_TISSUEMULTISCALEANALYSIS_H_
#include "bioLinearTissue.h"
#include "bioMultiscaleTissue.h"
#include "bioTissueAnalysis.h"
#include "bioVolumeConvergence.h"
#include <amsiMultiscale.h>
#include <amsiAnalysis.h>
#include <apfsimWrapper.h>
#include <apf.h>
#include <MeshSim.h>
#include <amsiPETScLAS.h>
#include <string>
#include <stdexcept>
namespace bio
{
  class MultiscaleTissueIteration : public amsi::Iteration
  {
  protected:
    MultiscaleTissue * tssu;
    amsi::LAS * las;
    Iteration * fem_iter;
  public:
  MultiscaleTissueIteration(MultiscaleTissue* a, amsi::LAS* l)
      : tssu(a), las(l), fem_iter(amsi::buildLinearFEMIteration(a, l))
  {
  }
  ~MultiscaleTissueIteration() { delete fem_iter; }
  void iterate();
  };
  class MultiscaleTissueAnalysis : public TissueAnalysis
  {
  public:
    MultiscaleTissueAnalysis(pGModel imdl, pParMesh imsh, pACase pd, MPI_Comm cm);
    virtual void init() override;
    virtual void run() override;
    virtual void finalizeStep() override;
  private:
    size_t cplng;
  };
} // end of namespace Biotissue
#endif
