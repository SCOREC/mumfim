#include "bioMultiscaleCoupling.h"
#include <apf.h>
#include "bioUtility.h"
#include <numeric>
namespace mumfim
{
  static double calcRVEDimensionality(apf::Mesh* mesh,
                               double fbr_area,
                               double fbr_vl_frc)
  {
    std::vector<double> lngths;
    calcDimMeasures(mesh,1,std::back_inserter(lngths));
    double ttl_fbr_lngth = std::accumulate(lngths.begin(),lngths.end(),0.0);
    return sqrt(ttl_fbr_lngth * fbr_area / fbr_vl_frc);
  }
  // FIXME the fiber area is the mean fiber area which isn't really a parameter (since
  // it is included in the mesh, and each fiber can have a different area)
  double calcScaleConversion(apf::Mesh* mesh, double fbr_area, double fbr_vol_frc)
  {
    double rve_dim = calcRVEDimensionality(mesh, fbr_area, fbr_vol_frc);
    double conversion_factor = 1/(rve_dim*rve_dim);
    return conversion_factor;
  }
}  // namespace mumfim
