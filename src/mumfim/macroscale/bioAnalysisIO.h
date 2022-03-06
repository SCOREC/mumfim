#ifndef BIO_ANALYSIS_IO_H_
#define BIO_ANALYSIS_IO_H_
#include <amsiUtil.h>
#include <apf.h>
namespace bio
{
  class NonlinearTissue;
  /**
   * @brief for each model item in the range [bgn_mdl_itm,nd_mdl_itm),
   *  print the step, entity tag, and volume to the specified log
   */
  template <typename I>
    void logVolumes(I bgn_mdl_itm, I nd_mdl_itm, amsi::Log log, int stp, apf::Field * U);
  /**
   * @brief for each model item in the range [bgn_mdl_itm,nd_mdl_tm),
   *  print the step, entity tag, and x, y, and z average displacement
   */
  template <typename I>
    void logDisps(I bgn_mdl_itm, I nd_mdl_itm, amsi::Log log, int stp, apf::Field * U);
  /**
   * @brief for each model item in the range [bgn, end), print the
   *  step, entity tag, and loads in the i, j, and k directions.
   */
  template <typename I>
    void logForces(I bgn, I end, amsi::Log log, int stp, NonlinearTissue * tssu);
}
#include "bioAnalysisIO_impl.h"
#endif
