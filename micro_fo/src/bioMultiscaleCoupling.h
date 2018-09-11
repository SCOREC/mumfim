#ifndef BIO_MULTISCALE_COUPLING_H_
#define BIO_MULTISCALE_COUPLING_H_
namespace apf
{
  class DynamicMatrix;
  class DynamicVector;
}  // namespace apf
namespace bio
{
  class FiberRVEAnalysis;
  class micro_fo_data;
  class micro_fo_result;
  class micro_fo_header;
  class micro_fo_params;
  class micro_fo_step_result;
  void recoverMultiscaleResults(FiberRVEAnalysis * ans, micro_fo_result * data);
  void recoverMultiscaleStepResults(FiberRVEAnalysis * ans,
                                    micro_fo_header & hdr,
                                    micro_fo_params & prm,
                                    micro_fo_step_result * data);
  void dsdxrve_2_dSdxrve(const apf::DynamicMatrix & ds_dx_rve,
                         const apf::DynamicVector & strs,
                         const apf::DynamicVector & dV_dx_rve,
                         double vol,
                         double cnv,
                         apf::DynamicMatrix & dS_dx_rve);
  void calcdR_dx_rve(apf::DynamicMatrix & dRdx_rve, FiberRVEAnalysis * ans);
  void calcdx_fn_dx_rve(apf::DynamicMatrix & dx_fn_dx_rve,
                        FiberRVEAnalysis * ans,
                        apf::DynamicMatrix & dR_dx_rve);
}  // namespace bio
#endif
