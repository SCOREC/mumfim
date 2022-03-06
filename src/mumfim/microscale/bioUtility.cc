#include "bioUtility.h"
namespace mumfim
{
  apf::Matrix3x3 eye()
  {
    apf::Matrix3x3 rslt;
    for (int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
        rslt[i][j] = (i==j) ? 1.0 : 0;
    return rslt;
  }
  apf::Matrix3x3 ones()
  {
    apf::Matrix3x3 rslt;
    for(int ii = 0; ii < 3; ii++)
      for(int jj = 0; jj < 3; jj++)
        rslt[ii][jj] = 1.0;
    return rslt;
  }
}
