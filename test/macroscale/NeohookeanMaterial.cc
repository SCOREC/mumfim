#include <catch2/catch.hpp>
#include "TestDeformations.h"
#include "mumfim/macroscale/materials/Materials.h"
#include "mumfim/macroscale/ApfMatrixUtils.h"
#include "TestSupport.h"
static apf::DynamicMatrix neohookean_stiffness(double shear_modulus,
                                               double poissons_ratio,
                                               apf::Matrix3x3 F)
{
  apf::DynamicMatrix material_stiffness(6,6);
  material_stiffness.zero();
  auto detF = apf::getDeterminant(F);
  double mu = shear_modulus;
  double lambda = (2.0 * mu * poissons_ratio) / (1.0 - 2.0 * poissons_ratio);
  auto lp = lambda/detF;
  auto mp = (mu-lambda * log(detF))/detF;
  material_stiffness(0,0) = lp + 2* mp;
  material_stiffness(1,1) = lp + 2* mp;
  material_stiffness(2,2) = lp + 2* mp;
  material_stiffness(3,3) = mp;
  material_stiffness(4,4) = mp;
  material_stiffness(5,5) = mp;
  material_stiffness(0,1) = material_stiffness(1,0) = lp;
  material_stiffness(0,2) = material_stiffness(2,0) = lp;
  material_stiffness(1,2) = material_stiffness(2,1) = lp;
  return material_stiffness;
}
static apf::DynamicMatrix neohookean_stress(double shear_modulus,
                                            double poissons_ratio,
                                            apf::Matrix3x3 F)
{
  double mu = shear_modulus;
  double lambda = (2.0 * mu * poissons_ratio) / (1.0 - 2.0 * poissons_ratio);
  apf::DynamicMatrix cauchy_stress(3,3);
  cauchy_stress.zero();
  auto b = mumfim::computeLeftCauchyGreen(F);
  auto detF = apf::getDeterminant(F);
  for(int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      cauchy_stress(i,j) = mu/detF*(b(i,j)-(i==j))+lambda/detF*log(detF)*(i==j);
    }
  }
  return cauchy_stress;
}
using test::compare_dynamic_matrices;
TEST_CASE("Neohookean Material", "[material][neohookean]")
{
  constexpr double poissons_ratio = 0.33;
  constexpr double shear_modulus = 100;
  double alpha = 0.5;
  SECTION("Pure Shear")
  {
    auto F = test::PureShear(0, alpha);
    const auto constitutive =
        mumfim::NeohookeanMaterial(shear_modulus, poissons_ratio, F);
    SECTION("Stress") {
    compare_dynamic_matrices(constitutive.cauchy_stress,
                             neohookean_stress(shear_modulus, poissons_ratio, F));
    }
    SECTION("Stiffness") {
      compare_dynamic_matrices(constitutive.material_stiffness,
                               neohookean_stiffness(shear_modulus, poissons_ratio, F));
    }
  }
  SECTION("Simple Shear")
  {
    auto F = test::SimpleShear(0, alpha);
    const auto constitutive =
        mumfim::NeohookeanMaterial(shear_modulus, poissons_ratio, F);
    SECTION("Stress") {
      compare_dynamic_matrices(constitutive.cauchy_stress,
                               neohookean_stress(shear_modulus, poissons_ratio, F));
    }
    SECTION("Stiffness") {
      compare_dynamic_matrices(constitutive.material_stiffness,
                               neohookean_stiffness(shear_modulus, poissons_ratio, F));
    }
  }
  SECTION("Uniaxial Extension")
  {
    auto F = test::UniaxialExtension(0, alpha);
    auto constitutive =
        mumfim::NeohookeanMaterial(shear_modulus, poissons_ratio, F);
    SECTION("Stress") {
      compare_dynamic_matrices(constitutive.cauchy_stress,
                               neohookean_stress(shear_modulus, poissons_ratio, F));
    }
    SECTION("Stiffness") {
      compare_dynamic_matrices(constitutive.material_stiffness,
                               neohookean_stiffness(shear_modulus, poissons_ratio, F));
    }
  }
}
