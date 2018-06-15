#ifndef BIO_MULTISCALE_RVE_ANALYSIS_IMPL_H_
#define BIO_MULTISCALE_RVE_ANALYSIS_IMPL_H_
#include "bioFiberRVEAnalysis.h"
#include <apfDynamicMatrix.h>
namespace bio
{
  void dsdxrve_2_dSdxrve(const apf::DynamicMatrix & ds_dx_rve,
                         const apf::DynamicVector & strs,
                         const apf::DynamicVector & dV_dx_rve,
                         double vol,
                         double cnv,
                         apf::DynamicMatrix & dS_dx_rve);
  void calcdR_dx_rve(apf::DynamicMatrix & dRdx_rve,
                     FiberRVEAnalysis * ans);
  void calcdx_fn_dx_rve(apf::DynamicMatrix & dx_fn_dx_rve,
                        FiberRVEAnalysis * ans,
                        apf::DynamicMatrix & dR_dx_rve);
}
#endif
