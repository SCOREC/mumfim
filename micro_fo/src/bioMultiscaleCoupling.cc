#include "bioMultiscaleCoupling.h"
#include <apf.h>
#include <apfMatrixUtil.h>
#include <cassert>
#include "bioFiberRVEAnalysis.h"
#include "bioMultiscaleMicroFOParams.h"
#include "bioUtil.h"
//#include "bioFiberNetwork.h"
#include <numeric>
namespace bio
{
  static double calcRVEDimensionality(FiberNetwork * fn,
                               double fbr_area,
                               double fbr_vl_frc)
  {
    std::vector<double> lngths;
    apf::Mesh * fn_msh = fn->getNetworkMesh();
    calcDimMeasures(fn_msh,1,std::back_inserter(lngths));
    double ttl_fbr_lngth = std::accumulate(lngths.begin(),lngths.end(),0.0);
    return sqrt(ttl_fbr_lngth * fbr_area / fbr_vl_frc);
  }
  double calcScaleConversion(FiberNetwork * fn, double fbr_area, double fbr_vol_frc)
  {
    double rve_dim = calcRVEDimensionality(fn, fbr_area, fbr_vol_frc);
    double conversion_factor = 1/(rve_dim*rve_dim);
    return conversion_factor;
  }
  void convertStressQuantities(FiberRVEAnalysis * ans, double * stress, double * C)
  {
    int dim = ans->getFn()->getNetworkMesh()->getDimension();
    int sigma_length = dim == 3 ? 6 : 3;
    int mat_stiff_length = dim == 3 ? 36 : 9;
    // convert to a macro-scale term
    double vol = ans->getRVE()->measureDu();
    double scale_conversion = ans->getFn()->getScaleConversion();
    for (int ii = 0; ii < sigma_length; ++ii)
      stress[ii] *= scale_conversion / vol;
    for (int ii = 0; ii < mat_stiff_length; ++ii)
      C[ii] *= scale_conversion / vol;
  }
  void recoverMultiscaleStepResults(RVEAnalysis * ans,
                                    micro_fo_header & hdr,
                                    micro_fo_params & prm,
                                    micro_fo_step_result * data)
  {
    FiberRVEAnalysis* FRveAns = dynamic_cast<FiberRVEAnalysis *>(ans);
    // if this if a fiber rve analysis set the orientation tensor
    if(FRveAns)
    {
      double * ornt_3d = &data->data[0];
      double * ornt_2d = &data->data[9];
      double n[3];
      n[0] = prm.data[ORIENTATION_AXIS_X];
      n[1] = prm.data[ORIENTATION_AXIS_Y];
      n[2] = prm.data[ORIENTATION_AXIS_Z];
      if (hdr.data[COMPUTE_ORIENTATION_3D])
      {
        get3DOrientationTensor(FRveAns->getFn(), ornt_3d);
      }
      if (hdr.data[COMPUTE_ORIENTATION_2D])
      {
        get2DOrientationTensor(FRveAns->getFn(), n, ornt_2d);
      }
    }
    // otherwise set the orientation tensor to 0 (since it has no meaning)
    else
    {
      for(int i=0; i<18; ++i)
        data->data[i] = 0;
    }
  }
}  // namespace bio
