#ifndef MUMFIM_TEST_MACROSCALE_TESTDEFORMATIONS_H
#define MUMFIM_TEST_MACROSCALE_TESTDEFORMATIONS_H
#include <apf.h>
#include <mumfim/exceptions.h>
namespace test
{
  [[nodiscard]] apf::Matrix3x3 Identity();
  [[nodiscard]] apf::Matrix3x3 PureShear(int direction = 0, double alpha = 1.0);
  [[nodiscard]] apf::Matrix3x3 SimpleShear(int direction = 0,
                                           double alpha = 1.0);
  [[nodiscard]] apf::Matrix3x3 UniaxialExtension(int direction = 0,
                                                 double alpha = 1.0);
}  // namespace test
#endif  // MUMFIM_TEST_MACROSCALE_TESTDEFORMATIONS_H
