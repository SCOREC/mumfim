#include "Materials.h"
#include "mumfim/macroscale/ApfMatrixUtils.h"

namespace mumfim {

  /*
    Determine stress-strain matrix, D, for NeoHookean material.
    See Bonet and Wood 2nd Edition, PP.250.
    - lambda and mu are effective lame parameters.
    - lambda = ShearModulus/detJ, mu = ( PoissonsRatio - ShearModulus *
    ln(detJ) )/detJ
  */
  // Calculate Cauchy stress from compressible NeoHookean equation. See
  // Bonet and Wood 2nd Ed. PP.163.
  MaterialResult NeohookeanMaterial(double shear_modulus,
                                    double poissons_ratio,
                                    const apf::Matrix3x3 & deformation_gradient)
  {
    double detF = getDeterminant(deformation_gradient);
    auto left_cauchy_green = computeLeftCauchyGreen(deformation_gradient);
    apf::DynamicMatrix cauchy_stress(3, 3);
    double mu = shear_modulus;
    double lambda = (2.0 * mu * poissons_ratio) / (1.0 - 2.0 * poissons_ratio);
    for (size_t i = 0; i < cauchy_stress.getRows(); i++)
    {
      for (size_t j = 0; j < cauchy_stress.getColumns(); j++)
      {
        cauchy_stress(i, j) = mu / detF * (left_cauchy_green(i, j) - (i == j)) +
            lambda / detF * log(detF) * (i == j);
      }
    }
    double lambda_prime = lambda / detF;
    double mu_prime = (mu - lambda * log(detF)) / detF;
    apf::DynamicMatrix D(6, 6);
    D.zero();
    D(0, 0) = lambda_prime + (2.0 * mu_prime);
    D(0, 1) = lambda_prime;
    D(0, 2) = lambda_prime;
    D(1, 0) = lambda_prime;
    D(1, 1) = lambda_prime + (2.0 * mu_prime);
    D(1, 2) = lambda_prime;
    D(2, 0) = lambda_prime;
    D(2, 1) = lambda_prime;
    D(2, 2) = lambda_prime + (2.0 * mu_prime);
    D(3, 3) = mu_prime;
    D(4, 4) = mu_prime;
    D(5, 5) = mu_prime;
    return MaterialResult{.cauchy_stress = cauchy_stress,
        .material_stiffness = D};
  }
}