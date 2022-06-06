#ifndef MUMFIM_SRC_MUMFIM_MACROSCALE_APFMATRIXUTILS_H
#define MUMFIM_SRC_MUMFIM_MACROSCALE_APFMATRIXUTILS_H
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>
namespace mumfim
{
  constexpr int num_symmetric_components(int n) { return (n * (n + 1)) / 2; }
  [[nodiscard]] apf::DynamicVector dynamic_matrix_to_voigt(
      const apf::DynamicMatrix & mat);
  [[nodiscard]] apf::Matrix3x3 computeGreenLagrangeStrain(
      const apf::Matrix3x3 & F);
  [[nodiscard]] apf::DynamicMatrix computeLeftCauchyGreen(
      const apf::Matrix3x3 & F);
  [[nodiscard]] apf::DynamicMatrix computeRightCauchyGreen(
      const apf::Matrix3x3 & F);
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_MACROSCALE_APFMATRIXUTILS_H
