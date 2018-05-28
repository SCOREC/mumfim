#ifndef BIO_MULTISCALE_RVE_ANALYSIS_IMPL_H_
#define BIO_MULTISCALE_RVE_ANALYSIS_IMPL_H_
#include "apfDynamicMatrix.h"
namespace bio
{
  void convertStressDiv(const apf::DynamicMatrix & ds_dx_rve,
                        const apf::DynamicVector & strs,
                        const apf::DynamicVector & dV_dx_rve,
                        double vol,
                        double cnv,
                        apf::DynamicMatrix & dS_dx_rve);
}
#endif
