#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_TEAM_OUTER_LOOP_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_TEAM_OUTER_LOOP_H__
#include "bioBatchedFiberRVEAnalysisExplicitBase.h"
namespace bio
{
  template <typename Scalar, typename Ordinal, typename ExeSpace>
  struct TeamOuterLoop
      : public BaseBatchedExplicit<TeamOuterLoop<Scalar, Ordinal, ExeSpace>,
                                   Scalar,
                                   Ordinal,
                                   ExeSpace>
  {
  };
}  // namespace bio
#endif
