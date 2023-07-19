#include <catch2/catch.hpp>

#include "TestDeformations.h"
#include "mumfim/exceptions.h"
#include "mumfim/macroscale/materials/Materials.h"

TEST_CASE("Load Torch", "[material][torch]")
{
  SECTION("Invalid file")
  {
    REQUIRE_THROWS_AS(mumfim::TorchMaterial("doesnt-exist.pt"),
                      mumfim::material_error);
  }
  SECTION("Valid file")
  {
    REQUIRE_NOTHROW(mumfim::TorchMaterial("data/kuhl-network.pt"));
  }
}

TEST_CASE("EvaluateTorch", "[material][torch]")
{
  auto material = mumfim::TorchMaterial("data/kuhl-network.pt");
  auto F = test::PureShear(0, 0.5);
  std::cout << F << std::endl;
  auto [cauchy_stress, stiffness] = material(F, nullptr, 0);
  std::cout << cauchy_stress << std::endl;
  std::cout << stiffness << std::endl;
  REQUIRE(cauchy_stress(0, 0) == Approx(75.0605));
  REQUIRE(cauchy_stress(1, 1) == Approx(75.0605));
  REQUIRE(cauchy_stress(2, 2) == Approx(18.9386));
  REQUIRE(cauchy_stress(1, 2) == Approx(0.0));
  REQUIRE(cauchy_stress(0, 2) == Approx(0.0));
  REQUIRE(cauchy_stress(0, 1) == Approx(79.8007));
  // VERIFY Cauchy stress is symmetric
  REQUIRE(cauchy_stress(1, 2) == Approx(cauchy_stress(2, 1)));
  REQUIRE(cauchy_stress(0, 2) == Approx(cauchy_stress(2, 0)));
  REQUIRE(cauchy_stress(0, 1) == Approx(cauchy_stress(1, 0)));

  // check the stiffness
  REQUIRE(stiffness(0, 0) == Approx(1.7782e+02));
  REQUIRE(stiffness(0, 1) == Approx(128.4382476807));
  REQUIRE(stiffness(0, 2) == Approx(2.8286e+01));
  REQUIRE(stiffness(0, 3) == Approx(-0.0000167617));
  REQUIRE(stiffness(0, 4) == Approx(-0.0000155636));
  REQUIRE(stiffness(0, 5) == Approx(1.4857e+02));

  // check stiffness is symmetric
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      // set a relatively large margin since pytorch is using single precision
      // the full stiffness tensor is constructed in pytorch (unlike the stress
      // which only fits the upper triangle) and then symmetrized in the pytorch
      // material
      REQUIRE(stiffness(i, j) == Approx(stiffness(j, i)).margin(1e-4));
    }
  }
  // tensor([[[ 1.7782e+02,  1.2844e+02,  2.8286e+01, -1.6766e-05,
  // -1.5566e-05, 1.4857e+02],
  //        [ 1.2844e+02,  1.7782e+02,  2.8286e+01, -1.5406e-05,
  //        -1.6206e-05, 1.4857e+02], [ 2.8286e+01,  2.8286e+01,  8.8807e+01,
  //        -1.0790e-05, -1.0790e-05, 3.5535e+01], [ 0.0000e+00,  0.0000e+00,
  //        0.0000e+00,  5.1216e+01,  3.8583e+01, 0.0000e+00], [ 0.0000e+00,
  //        0.0000e+00,  0.0000e+00,  3.8583e+01,  5.1216e+01, 0.0000e+00],
  //        [ 1.4857e+02,  1.4857e+02,  3.5535e+01, -1.2869e-05,
  //        -1.2709e-05, 1.4618e+02]]],
}