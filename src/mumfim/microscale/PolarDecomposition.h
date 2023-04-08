#ifndef POLARDECOMP_POLAR_DECOMPOSITION_H
#define POLARDECOMP_POLAR_DECOMPOSITION_H
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Gemv_Serial_Impl.hpp>
#include <KokkosBatched_SVD_Decl.hpp>
#include <KokkosBatched_Set_Decl.hpp>
#include <KokkosBatched_Set_Impl.hpp>
#include <Kokkos_Core.hpp>
namespace mumfim
{
  template <int N>
  struct PolarResult
  {
    PolarResult(int size) : P("P", size), U("U", size) {}
    /**
     * Hermetian positive semidefinite. If a is nonsingular, P is positive
     * definite
     */
    Kokkos::View<double * [N][N]> P;
    /**
     * Unitary matrix
     */
    Kokkos::View<double * [N][N]> U;
  };
  template <>
  struct PolarResult<0>
  {
    PolarResult(int size, int N) : P("P", size, N, N), U("U", size, N, N) {}
    Kokkos::View<double ***> P;
    Kokkos::View<double ***> U;
  };
  struct LeftTag
  {
  };
  struct RightTag
  {
  };
  /**
   * @tparam Side Determines side (left/right of polar decomposition If
   * side=LeftTag then a=up. If side=RightTag a=pu.
   * @tparam N size of the matrix
   * @param matrix the matrix that the polar decomposition should be computed on
   * @return PolarResult
   */
  template <typename ExeSpace, typename Side = RightTag, int N>
  auto PolarDecomposition(
      Kokkos::View<double * [N][N], typename ExeSpace::memory_space> matrices)
      -> PolarResult<N>
  {
    static_assert(std::is_same_v<Side, RightTag>,
                  "Currently only right polar decomposition implemented");
    Kokkos::View<double * [N][N]> Us("Us", matrices.extent(0));
    Kokkos::View<double * [N]> Ss("Ss", matrices.extent(0));
    Kokkos::View<double * [N][N]> Vts("Vts", matrices.extent(0));
    Kokkos::View<double * [N]> works("works", matrices.extent(0));
    Kokkos::View<double * [N][N]> matrices_copy("matrices_copy",
                                                matrices.extent(0));
    Kokkos::deep_copy(matrices_copy, matrices);
    PolarResult<N> results(matrices.extent(0));
    auto policy = Kokkos::RangePolicy<ExeSpace>(0, matrices.extent(0));
    Kokkos::parallel_for(
        "polar decomposition", policy, KOKKOS_LAMBDA(int i) {
          auto matrix_copy =
              Kokkos::subview(matrices_copy, i, Kokkos::ALL(), Kokkos::ALL());
          auto U = Kokkos::subview(Us, i, Kokkos::ALL(), Kokkos::ALL());
          auto S = Kokkos::subview(Ss, i, Kokkos::ALL());
          auto Vt = Kokkos::subview(Vts, i, Kokkos::ALL(), Kokkos::ALL());
          auto work = Kokkos::subview(works, i, Kokkos::ALL());
          auto result_U =
              Kokkos::subview(results.U, i, Kokkos::ALL(), Kokkos::ALL());
          auto result_P =
              Kokkos::subview(results.P, i, Kokkos::ALL(), Kokkos::ALL());
          KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag{},
                                           matrix_copy, U, S, Vt, work);
          KokkosBatched::SerialGemm<
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, U, Vt, 0.0,
                                                            result_U);
          for (size_t k = 0; k < U.extent(0); ++k)
          {
            for (size_t l = 0; l < U.extent(1); ++l)
            {
              U(k, l) = Vt(l, k) * S(l);
            }
          }
          KokkosBatched::SerialGemm<
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, U, Vt, 0.0,
                                                            result_P);
        });
    return results;
  }
}  // namespace mumfim
#endif  // POLARDECOMP_POLAR_DECOMPOSITION_H
