#ifndef BIO_MULTISCALE_COUPLING_H_
#define BIO_MULTISCALE_COUPLING_H_
namespace apf
{
  class DynamicMatrix;
  class DynamicVector;
  class Mesh;
}  // namespace apf
namespace bio
{
  class RVEAnalysis;
  class micro_fo_data;
  class micro_fo_result;
  class micro_fo_header;
  class micro_fo_params;
  class micro_fo_step_result;
  class FiberNetwork;
  class DeformationGradient;
  class FiberRVEAnalysis;
  void recoverMultiscaleStepResults(RVEAnalysis * ans,
                                    micro_fo_header & hdr,
                                    micro_fo_params & prm,
                                    micro_fo_step_result * data);
  /*
   * \brief computes RVE scaling factor
   * Computes the scaling factor. This correlates to the RVE side length in "physical space"
   * \param fn a pointer to the fiber network data structure
   * \param fbr_area fiber cross-sectional area (from experiment)
   * \fbr_vol_frc volume fraction of the fiber network (from experiment)
   * \warning This function assumes that all of the fibers have the same cross-sectional area
   */
  double calcScaleConversion(apf::Mesh * mesh, double fbr_area, double fbr_vol_frc);
  void convertStressQuantities(FiberRVEAnalysis * ans, double * stress, double * C);
}  // namespace bio
#endif
