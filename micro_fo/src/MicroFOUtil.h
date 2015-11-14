#ifndef BIO_MICRO_FO_UTIL_H_
#define BIO_MICRO_FO_UTIL_H_

#include <cmath>
#include <limits>

namespace bio
{
  inline bool close(double x, double y)
  {
    return fabs(x - y) < std::numeric_limits<double>::epsilon();
  }
}
#endif
