#ifndef MUMFIM_SINGLESCALETISSUEANALYSIS_H
#define MUMFIM_SINGLESCALETISSUEANALYSIS_H
#include "TissueAnalysis.h"

namespace mumfim
{
  class SinglescaleTissueAnalysis : public TissueAnalysis
  {
    public:
    SinglescaleTissueAnalysis(apf::Mesh * mesh,
    std::unique_ptr<const mt::CategoryNode> cs,
        MPI_Comm c,
    const amsi::Analysis & amsi_analysis);
  };
}  // namespace mumfim
#endif
