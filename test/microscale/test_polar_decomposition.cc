#include <mumfim/microscale/PolarDecomposition.h>
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACKE
#include <KokkosBatched_Eigendecomposition_Decl.hpp>
#include <KokkosBatched_Eigendecomposition_Serial_Impl.hpp>
#endif
#include <catch2/catch.hpp>
#include <KokkosBlas3_gemm.hpp>


TEST_CASE("Polar Decomposition") {
  using exe_space = Kokkos::Serial;
  using memory_space = typename exe_space::memory_space;
  constexpr int num_trials = 5;
  Kokkos::View<double*[3][3], memory_space> matrix("mat",num_trials);
  Kokkos::deep_copy(matrix, 0);
  Kokkos::parallel_for("fill",Kokkos::RangePolicy<exe_space>(0, num_trials) , KOKKOS_LAMBDA(int i) {
      matrix(i, 0, 0) = 1;
      matrix(i, 0, 1) = 8;
      matrix(i, 1, 0) = 0.2;
      matrix(i, 1, 1) = 0.7;
      matrix(i, 2, 2) = 0.4;
      });
  auto [P, U] = mumfim::PolarDecomposition<exe_space>(matrix);
  //Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, 
  // verify rank of P and U is 2
  REQUIRE(P.rank == 3);
  REQUIRE(U.rank == 3);
  // verify extents of P and U are 3x3
  REQUIRE(P.extent(0) == num_trials);
  REQUIRE(P.extent(1) == 3);
  REQUIRE(P.extent(2) == 3);
  REQUIRE(U.extent(0) == num_trials);
  REQUIRE(U.extent(1) == 3);
  REQUIRE(U.extent(2) == 3);


  auto P_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, P);
  auto U_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, U);
  auto P0 = Kokkos::subview(P_h, 1, Kokkos::ALL(), Kokkos::ALL());
  auto U0 = Kokkos::subview(U_h, 1, Kokkos::ALL(), Kokkos::ALL());
  // verify U is unitary
  SECTION("U is unitary") {
    Kokkos::View<double[3][3], memory_space> result("result");
    KokkosBlas::gemm("N", "T", 1.0, U0, U0,0, result);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        std::cout << U0(i, j) << " ";
        if (i == j) {
          REQUIRE(result(i, j) == Approx(1.0));
        } else {
          REQUIRE(result(i, j) == Approx(0.0));
        }
      }
      std::cout<<"\n";
    }
  }
  SECTION("P is symmetric") {
    // verify that P is symmetric
    for (size_t i = 0; i < P0.extent(0); ++i) {
      for (size_t j = 0; j < P0.extent(1); ++j) {
        std::cout << P0(i, j) << " ";
        REQUIRE(P0(i, j) == Approx(P0(j, i)));
      }
      std::cout <<"\n";
    }
  }
  SECTION("UP=matrix") {
    auto matrix_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrix);
    // verify that U*P = matrix
    Kokkos::View<double[3][3],memory_space> result("result");
    KokkosBlas::gemm("N", "N", 1.0, U0, P0, 0, result);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        REQUIRE(result(i, j) == Approx(matrix_h(0,i, j)));
      }
    }
  }
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACKE
  SECTION("P is positive definite") {
    // verify that P is positive definite
    Kokkos::View<double *,memory_space> er("er", P.extent(0));
    Kokkos::View<double *,memory_space> ei("ei", P.extent(0));
    Kokkos::View<double **, memory_space> UL("UL", P.extent(0), P.extent(1));
    Kokkos::View<double **, memory_space> UR("UR", P.extent(0), P.extent(1));
    Kokkos::View<double *, memory_space> work("work", 2 * P.extent(0) * P.extent(0) + 5 * P.extent(0));
    KokkosBatched::SerialEigendecomposition::invoke(P, er, ei, UL, UR, work);
    for (int i = 0; i < P.extent(0); ++i) {
      REQUIRE(er(i) > 0);
    }
  }
#endif
}
