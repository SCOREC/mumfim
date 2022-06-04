#include "SinglescaleTissueAnalysis.h"
mumfim::SinglescaleTissueAnalysis::SinglescaleTissueAnalysis(
    apf::Mesh * mesh,
    std::unique_ptr<const mt::CategoryNode> cs,
    MPI_Comm c,
    const amsi::Analysis & amsi_analysis)
    : TissueAnalysis(mesh, std::move(cs), c, amsi_analysis)
{
  const auto * solution_strategy =
      mt::GetPrimaryCategoryByType(analysis_case.get(), "solution strategy");
  tssu = new NonlinearTissue(mesh, *analysis_case, cm);
  addVolumeTracking(mesh, solution_strategy);
  // We want to do the tissue iteration after we compute the volumes
  itr_stps.push_back(new TissueIteration(tssu, las));
  itr_stps.push_back(new TissueCheckpointIteration(this));
  itr = new amsi::MultiIteration(itr_stps.begin(), itr_stps.end());
  buildLASConvergenceOperators(solution_strategy,itr,las,std::back_inserter(cvg_stps));
  buildVolConvergenceOperators(solution_strategy,itr,tssu->getUField(),trkd_vols,std::back_inserter(cvg_stps));
  cvg = new amsi::MultiConvergence(cvg_stps.begin(), cvg_stps.end());
}
