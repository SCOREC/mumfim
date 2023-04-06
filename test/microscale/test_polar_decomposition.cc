#include <mumfim/microscale/PolarDecomposition.h>
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACKE
#include <KokkosBatched_Eigendecomposition_Decl.hpp>
#include <KokkosBatched_Eigendecomposition_Serial_Impl.hpp>
#endif
//#include <Kokkos_core.hpp>
#include <catch2/catch.hpp>

TEST_CASE("Polar Decomposition") {
  Kokkos::View<double[3][3]> matrix("mat");
  Kokkos::deep_copy(matrix, 0);
  matrix(0, 0) = 1;
  matrix(0, 1) = 8;
  matrix(1, 0) = 0.2;
  matrix(1, 1) = 0.7;
  matrix(2, 2) = 0.4;
  auto [P, U] = polar_decomposition(matrix);
  // polar decomposition should not modify the input matrix
  REQUIRE(matrix(0,0) == Approx(1.0));
  REQUIRE(matrix(0,1) == Approx(8.0));
  REQUIRE(matrix(1,0) == Approx(0.2));
  REQUIRE(matrix(1,1) == Approx(0.7));
  REQUIRE(matrix(2,2) == Approx(0.4));
  // verify rank of P and U is 2
  REQUIRE(P.rank == 2);
  REQUIRE(U.rank == 2);
  // verify extents of P and U are 3x3
  REQUIRE(P.extent(0) == 3);
  REQUIRE(P.extent(1) == 3);
  REQUIRE(U.extent(0) == 3);
  REQUIRE(U.extent(1) == 3);
  // verify U is unitary
  SECTION("U is unitary") {
    Kokkos::View<double[3][3]> result("result");
    KokkosBlas::gemm("N", "T", 1.0, U,U,0, result);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (i == j) {
          REQUIRE(result(i, j) == Approx(1.0));
        } else {
          REQUIRE(result(i, j) == Approx(0.0));
        }
      }
    }
  }
  SECTION("P is symmetric") {
    // verify that P is symmetric
    for (size_t i = 0; i < P.extent(0); ++i) {
      for (size_t j = 0; j < P.extent(1); ++j) {
        REQUIRE(P(i, j) == Approx(P(j, i)));
      }
    }
  }
  SECTION("UP=matrix") {
    // verify that U*P = matrix
    Kokkos::View<double[3][3]> result("result");
    KokkosBlas::gemm("N", "N", 1.0, U, P, 0, result);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        REQUIRE(result(i, j) == Approx(matrix(i, j)));
      }
    }
  }
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACKE
  SECTION("P is positive definite") {
    // verify that P is positive definite
    Kokkos::View<double *> er("er", P.extent(0));
    Kokkos::View<double *> ei("ei", P.extent(0));
    Kokkos::View<double **> UL("UL", P.extent(0), P.extent(1));
    Kokkos::View<double **> UR("UR", P.extent(0), P.extent(1));
    Kokkos::View<double *> work("work", 2 * P.extent(0) * P.extent(0) + 5 * P.extent(0));
    KokkosBatched::SerialEigendecomposition::invoke(P, er, ei, UL, UR, work);
    for (int i = 0; i < P.extent(0); ++i) {
      REQUIRE(er(i) > 0);
    }
  }
#endif
}