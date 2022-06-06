#include "TestDeformations.h"
namespace test
{
  apf::Matrix3x3 Identity()
  {
    apf::Matrix3x3 F;
    for(int i=0; i<3; ++i) {
      for(int j=0; j<3; ++j){
        F[i][j] = 0;
      }
    }
    for (int i = 0; i < 3; ++i)
    {
      F[i][i] = 1;
    }
    return F;
  }
  apf::Matrix3x3 PureShear(int direction, double alpha)
  {
    if (direction < 0 || direction > 2)
    {
      throw mumfim::mumfim_error("direction must be between 0 and 2");
    }
    auto F = Identity();
    if (direction == 0)
    {
      F[0][1] = F[1][0] = alpha;
    }
    else if (direction == 1)
    {
      F[0][2] = F[2][0] = alpha;
    }
    else if (direction == 2)
    {
      F[1][2] = F[2][1] = alpha;
    }
    return F;
  }
  apf::Matrix3x3 SimpleShear(int direction, double alpha)
  {
    if (direction < 0 || direction > 5)
    {
      throw mumfim::mumfim_error("direction must be between 0 and 6");
    }
    auto F = Identity();
    if (direction == 0)
    {
      F[0][1] = alpha;
    }
    else if (direction == 1)
    {
      F[0][2] = alpha;
    }
    else if (direction == 2)
    {
      F[1][2] = alpha;
    }
    if (direction == 3)
    {
      F[1][0] = alpha;
    }
    else if (direction == 4)
    {
      F[2][0] = alpha;
    }
    else if (direction == 5)
    {
      F[2][1] = alpha;
    }
    return F;
  }
  apf::Matrix3x3 UniaxialExtension(int direction, double alpha)
  {
    if (direction < 0 || direction > 2)
    {
      throw mumfim::mumfim_error("direction must be between 0 and 2");
    }
    auto F = Identity();
    F[direction][direction] += alpha;
    return F;
  }
}  // namespace test
