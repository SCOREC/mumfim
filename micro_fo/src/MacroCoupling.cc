#include "MacroCoupling.h"
#include "FiberNetwork.h"

#include <numeric>
#include <vector>

namespace bio
{
  double calcRVEDimensionality(const FiberNetwork * fn,double fbr_area, double fbr_vl_frc)
  {
    std::vector<double> fbr_lngths;
    fn->calcFiberLengths(fbr_lngths);
    double ttl_fbr_lngth = std::accumulate(fbr_lngths.begin(),
					  fbr_lngths.end(),
					  0.0);
    return sqrt(ttl_fbr_lngth * fbr_area / fbr_vl_frc);
  }
}
