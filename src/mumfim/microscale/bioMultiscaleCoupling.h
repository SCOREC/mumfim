#ifndef MUMFIM_MULTISCALE_COUPLING_H_
#define MUMFIM_MULTISCALE_COUPLING_H_
#include <apf.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include "bioBatchedFiberRVEAnalysisExplicit.h"
#include "bioMultiscaleMicroFOParams.h"
namespace mumfim
{
  template <typename DualView>
  void copyOrientationResultsFromKokkosArray(
      DualView orientation_tensor,
      std::vector<micro_fo_step_result> & results,
      int offset)
  {
    orientation_tensor.template sync<Kokkos::HostSpace>();
    auto orientation_tensor_h = orientation_tensor.h_view;
    for (std::size_t i = 0; i < results.size(); ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        for (int k = 0; k < 3; ++k)
        {
          results[i].data[j * 3 + k + offset] = orientation_tensor_h(i, j, k);
        }
      }
    }
  }
  template <typename DualView,
            typename DualViewNorm,
            typename BatchedAnalysisType>
  void recoverMultiscaleStepResults(DualView orientation_tensor,
                                    DualViewNorm orientation_normal,
                                    BatchedAnalysisType * batched_analysis,
                                    std::vector<micro_fo_header> & /*hdrs*/,
                                    std::vector<micro_fo_params> & prms,
                                    std::vector<micro_fo_step_result> & results)
  {
    orientation_normal.template modify<Kokkos::HostSpace>();
    for (std::size_t i = 0; i < results.size(); ++i)
    {
      orientation_normal.h_view(i, 0) = prms[i].data[ORIENTATION_AXIS_X];
      orientation_normal.h_view(i, 1) = prms[i].data[ORIENTATION_AXIS_Y];
      orientation_normal.h_view(i, 2) = prms[i].data[ORIENTATION_AXIS_Z];
    }
    batched_analysis->compute3DOrientationTensor(orientation_tensor);
    copyOrientationResultsFromKokkosArray(orientation_tensor, results, 0);
    batched_analysis->compute2DOrientationTensor(orientation_normal,
                                                 orientation_tensor);
    copyOrientationResultsFromKokkosArray(orientation_tensor, results, 9);
  }
  /*
   * \brief computes RVE scaling factor
   * Computes the scaling factor. This correlates to the RVE side length in "physical space"
   * \param fn a pointer to the fiber network data structure
   * \param fbr_area fiber cross-sectional area (from experiment)
   * \fbr_vol_frc volume fraction of the fiber network (from experiment)
   * \warning This function assumes that all of the fibers have the same cross-sectional area
   */
  double calcScaleConversion(apf::Mesh * mesh, double fbr_area, double fbr_vol_frc);
}  // namespace mumfim
#endif
