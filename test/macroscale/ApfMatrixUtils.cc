#include <catch2/catch.hpp>
#include <mumfim/macroscale/ApfMatrixUtils.h>
#include "TestDeformations.h"

TEST_CASE("Grean Lagrange Strain") {
  SECTION("Pure Shear") {
    double alpha = 2.0;
    for (int i=0;i<4;++i)
    {
      const auto F = test::PureShear(0, alpha);
      const auto GL0 = mumfim::computeGreenLagrangeStrain(F);
      REQUIRE(GL0[0][0] == Approx(1/2.0*alpha*alpha));
      REQUIRE(GL0[0][1] == Approx( alpha));
      REQUIRE(GL0[0][2] == Approx(0));
      REQUIRE(GL0[1][0] == Approx( alpha));
      REQUIRE(GL0[1][1] == Approx(1/2.0*alpha*alpha));
      REQUIRE(GL0[1][2] == Approx(0));
      REQUIRE(GL0[2][0] == Approx(0));
      REQUIRE(GL0[2][1] == Approx(0));
      REQUIRE(GL0[2][2] == Approx(0));
      alpha *= 2.0;
    }
  }
  SECTION("Simple Shear") {
    double alpha = 2.0;
    for (int i=0;i<4;++i) {
      const auto F = test::SimpleShear(0,alpha);
      const auto GL0 = mumfim::computeGreenLagrangeStrain(F);
      REQUIRE(GL0[0][0] == Approx(0));
      REQUIRE(GL0[0][1] == Approx(alpha/2.0));
      REQUIRE(GL0[0][2] == Approx(0));
      REQUIRE(GL0[1][0] == Approx(alpha/2.0));
      REQUIRE(GL0[1][1] == Approx(alpha*alpha/2.0));
      REQUIRE(GL0[1][2] == Approx(0));
      REQUIRE(GL0[2][0] == Approx(0));
      REQUIRE(GL0[2][1] == Approx(0));
      REQUIRE(GL0[2][2] == Approx(0));
      alpha *= 2.0;
    }
  }
  SECTION("Uniaxial Extension") {
    double alpha = 2.0;
    for (int i=0;i<4;++i) {
      const auto F = test::UniaxialExtension(0,alpha);
      const auto GL0 = mumfim::computeGreenLagrangeStrain(F);
      REQUIRE(GL0[0][0] == Approx(((1+alpha)*(1+alpha)-1)/2.0));
      REQUIRE(GL0[0][1] == Approx(0));
      REQUIRE(GL0[0][2] == Approx(0));
      REQUIRE(GL0[1][0] == Approx(0));
      REQUIRE(GL0[1][1] == Approx(0));
      REQUIRE(GL0[1][2] == Approx(0));
      REQUIRE(GL0[2][0] == Approx(0));
      REQUIRE(GL0[2][1] == Approx(0));
      REQUIRE(GL0[2][2] == Approx(0));
      alpha *= 2.0;
    }
  }
}