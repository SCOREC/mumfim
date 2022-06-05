#include "UpdatedLagrangianMaterialIntegrator.h"
#include <ElementalSystem.h>

namespace mumfim
{
  apf::DynamicVector dynamic_matrix_to_voigt(const apf::DynamicMatrix & mat)
  {
    if (mat.getRows() != mat.getColumns() || mat.getRows() != 3 ||
        mat.getRows() != 2)
    {
      throw mumfim_error("mat should be 3x3 or 2x2 matrix");
    }
    apf::DynamicVector vec(num_symmetric_components(mat.getColumns()));
    for (int i; i < mat.getColumns(); ++i)
    {
      vec(i) = mat(i, i);
    }
    if (mat.getColumns() == 3)
    {
      vec(3) = mat(2, 3);
      vec(4) = mat(1, 3);
      vec(5) = mat(1, 2);
    }
    else if (mat.getColumns() == 2)
    {
      vec(3) = mat(1, 2);
    }
  }
  apf::Matrix3x3 computeGreenLagrangeStrain(const apf::Matrix3x3 & F)
  {
    auto rightCauchyGreen = computeRightCauchyGreen(F);
    // E_G = 1/2(C-I), C=F^T.F, Green-Lagrange Strain.
    return {0.5 * (rightCauchyGreen(0, 0) - 1), 0.5 * rightCauchyGreen(0, 1),
            0.5 * rightCauchyGreen(0, 2),       0.5 * rightCauchyGreen(1, 0),
            0.5 * (rightCauchyGreen(1, 1) - 1), 0.5 * rightCauchyGreen(1, 2),
            0.5 * rightCauchyGreen(2, 0),       0.5 * rightCauchyGreen(2, 1),
            0.5 * (rightCauchyGreen(2, 2) - 1)};
  }
  apf::DynamicMatrix computeLeftCauchyGreen(const apf::Matrix3x3 & F) {

    apf::DynamicMatrix leftCauchyGreen(3, 3);  // rightCauchyGreen.zero();
    apf::DynamicMatrix FT(3, 3);
    FT.zero();
    apf::transpose(fromMatrix(F), FT);
    apf::multiply(fromMatrix(F),FT, leftCauchyGreen);
    return leftCauchyGreen;
  }
  apf::DynamicMatrix computeRightCauchyGreen(const apf::Matrix3x3 & F) {

    apf::DynamicMatrix rightCauchyGreen(3, 3);  // rightCauchyGreen.zero();
    apf::DynamicMatrix FT(3, 3);
    FT.zero();
    apf::transpose(fromMatrix(F), FT);
    apf::multiply(FT, fromMatrix(F), rightCauchyGreen);
    return rightCauchyGreen;
  }
}  // namespace mumfim
