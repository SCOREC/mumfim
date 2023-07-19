#ifndef POLARDECOMP_POLAR_DECOMPOSITION_H
#define POLARDECOMP_POLAR_DECOMPOSITION_H
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Gemv_Serial_Impl.hpp>
#include <KokkosBatched_SVD_Decl.hpp>
#include <KokkosBatched_Set_Decl.hpp>
// #include <KokkosBatched_Set_Impl.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Copy_Impl.hpp>

namespace mumfim
{
  template <int N, typename MemorySpace>
  struct PolarResult
  {
    PolarResult(int size) : U("P", size), R("U", size) {}
    /**
     * Hermetian positive semidefinite. If a is nonsingular, U is positive
     * definite (also often called P)
     */
    Kokkos::View<double * [N][N], MemorySpace> U;
    /**
     * Unitary matrix (also often called U)
     */
    Kokkos::View<double * [N][N], MemorySpace> R;
  };
  template <typename MemorySpace>
  struct PolarResult<0,MemorySpace>
  {
    PolarResult(int size, int N) : U("U", size, N, N), R("R", size, N, N) {}
    Kokkos::View<double ***,MemorySpace> U;
    Kokkos::View<double ***,MemorySpace> R;
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
  //template <typename ExeSpace, typename Side = RightTag, int N>
  template <typename ExeSpace, typename ViewT, int N=3, typename Side=RightTag>
  auto PolarDecomposition(
      //Kokkos::View<double * [N][N], typename ExeSpace::memory_space> matrices)
      ViewT matrices)
      -> PolarResult<N, typename ExeSpace::memory_space>
  {
    using memory_space = typename ExeSpace::memory_space;
    static_assert(std::is_same_v<Side, RightTag>,
                  "Currently only right polar decomposition implemented");
    // SVD work arrays need to be layout right or the algorithm is broken see https://github.com/kokkos/kokkos-kernels/issues/1786
    using array_layout = Kokkos::LayoutRight;
    constexpr bool layouts_match = std::is_same_v<array_layout,typename ViewT::array_layout>;
    Kokkos::View<double * [N][N], array_layout, memory_space> Us("Us", matrices.extent(0));
    Kokkos::View<double * [N], array_layout, memory_space> Ss("Ss", matrices.extent(0));
    Kokkos::View<double * [N][N],array_layout, memory_space> Vts("Vts", matrices.extent(0));
    Kokkos::View<double * [N],array_layout, memory_space> works("works", matrices.extent(0));
    Kokkos::View<double * [N][N], array_layout, memory_space> matrices_copy("matrices_copy",
                                                matrices.extent(0));
    PolarResult<N,typename ExeSpace::memory_space> results(matrices.extent(0));
    Kokkos::View<double * [N][N], array_layout, memory_space> U_result;
    Kokkos::View<double * [N][N], array_layout, memory_space> R_result;
    if constexpr (layouts_match) {
      U_result = results.U;
      R_result = results.R;
    }
    else {
      U_result = Kokkos::View<double * [N][N], array_layout, memory_space>("U", matrices.extent(0));
      R_result = Kokkos::View<double * [N][N], array_layout, memory_space>("R", matrices.extent(0));
    }
    Kokkos::deep_copy(matrices_copy, matrices);
    auto policy = Kokkos::RangePolicy<ExeSpace>(0, matrices.extent(0));
    Kokkos::parallel_for(
        "polar decomposition", policy, KOKKOS_LAMBDA(int i) {
          auto matrix_copy =
              Kokkos::subview(matrices_copy, i, Kokkos::ALL(), Kokkos::ALL());
          auto U = Kokkos::subview(Us, i, Kokkos::ALL(), Kokkos::ALL());
          auto S = Kokkos::subview(Ss, i, Kokkos::ALL());
          auto Vt = Kokkos::subview(Vts, i, Kokkos::ALL(), Kokkos::ALL());
          auto work = Kokkos::subview(works, i, Kokkos::ALL());
          auto result_R =
              Kokkos::subview(R_result, i, Kokkos::ALL(), Kokkos::ALL());
          auto result_U =
              Kokkos::subview(U_result, i, Kokkos::ALL(), Kokkos::ALL());
          KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag{},
                                           matrix_copy, U, S, Vt, work);
          KokkosBatched::SerialGemm<
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, U, Vt, 0.0,
                                                            result_R);
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
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, U, Vt, 0.0, result_U);
        });
    if constexpr (! layouts_match) {
      Kokkos::deep_copy(results.R, R_result);
      Kokkos::deep_copy(results.U, U_result);
    }
    return results;
  }
}  // namespace mumfim
#endif  // POLARDECOMP_POLAR_DECOMPOSITION_H
