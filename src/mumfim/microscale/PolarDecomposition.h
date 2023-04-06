#ifndef POLARDECOMP_POLAR_DECOMPOSITION_H
#define POLARDECOMP_POLAR_DECOMPOSITION_H
#include <KokkosBatched_SVD_Decl.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_Core.hpp>

template<int N>
struct PolarResult {
  PolarResult() : P("P"), U("U") {}
  /**
   * Hermetian positive semidefinite. If a is nonsingular, P is positive definite
   */
  Kokkos::View<double[N][N]> P;
  /**
   * Unitary matrix
   */
  Kokkos::View<double[N][N]> U;
};

template<>
struct PolarResult<0> {
  PolarResult(int N) : P("P", N, N), U("U", N, N) {}
  Kokkos::View<double **> P;
  Kokkos::View<double **> U;
};

struct LeftTag {};
struct RightTag {};

/**
 * @tparam Side Determines side (left/right of polar decomposition If side=LeftTag then a=up. If side=RightTag a=pu.
 * @tparam N size of the matrix
 * @param matrix the matrix that the polar decomposition should be computed on
 * @return PolarResult
 */
template<typename Side = RightTag, int N>
PolarResult<N> polar_decomposition(Kokkos::View<double[N][N]> matrix) {
  static_assert(std::is_same_v<Side, RightTag>, "Currently only right polar decomposition implemented");
  Kokkos::View<double[N][N]> U("U");
  Kokkos::View<double[N]> S("S");
  Kokkos::View<double[N][N]> Vt("Vt");
  Kokkos::View<double[N]> work("work");
  Kokkos::View<double[N][N]> matrix_copy("matrix_copy");
  Kokkos::deep_copy(matrix_copy, matrix);
  KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag{}, matrix_copy, U, S, Vt, work);
  PolarResult<N> result;
  KokkosBlas::gemm("N", "N", 1.0, U, Vt, 0, result.U);
  Kokkos::View<double[N][N]> v;
  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>> {{0,0},{N,N} }, KOKKOS_LAMBDA(int i, int j) {
        U(i, j) = Vt(j, i) * S(j);
      });
  KokkosBlas::gemm("N", "N", 1.0, U, Vt, 0, result.P);
  return result;
}

#endif//POLARDECOMP_POLAR_DECOMPOSITION_H
