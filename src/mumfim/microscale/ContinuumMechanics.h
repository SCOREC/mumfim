#ifndef MUMFIM_SRC_MUMFIM_MICROSCALE_STRESSCONVERSION_H
#define MUMFIM_SRC_MUMFIM_MICROSCALE_STRESSCONVERSION_H
#include <Kokkos_Core.hpp>
#include "mumfim/microscale/MicroTypeDefinitions.h"
#include "mumfim/microscale/TensorUtilities.h"
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_SetIdentity_Decl.hpp"
#include "KokkosBatched_SetIdentity_Impl.hpp"

namespace mumfim {
  /**
   * @brief convert the stress tensor from Cauchy to PK2
   * @param[in] F the deformation gradient
   * @param[in,out] stress the Cauchy stress tensor (in Voigt form)
   */
  template <typename ExeSpace>
  void ConvertCauchyToPK2(
      Kokkos::View<Scalar * [3][3], typename ExeSpace::memory_space> F,
      Kokkos::View<Scalar * [6], typename ExeSpace::memory_space> stress)
  {
    Kokkos::View<Scalar * [3][3], typename ExeSpace::memory_space> invF(
        "invF", F.extent(0));
    Kokkos::View<Scalar * [3][3], typename ExeSpace::memory_space>
        stress_matrix("stress matrix", F.extent(0));
    // a work array
    Kokkos::View<Scalar * [3][3], typename ExeSpace::memory_space> work(
        "work", F.extent(0));
    assert(F.extent(0) == stress.extent(0));
    Kokkos::parallel_for(
        "convert cauchy to pk2",
        Kokkos::RangePolicy<ExeSpace>(0, stress.extent(0)),
        KOKKOS_LAMBDA(int i) {
      auto invFi = Kokkos::subview(invF, i, Kokkos::ALL(), Kokkos::ALL());
      auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
      auto StressVoigti = Kokkos::subview(stress, i, Kokkos::ALL());
      auto Stressi =
          Kokkos::subview(stress_matrix, i, Kokkos::ALL(), Kokkos::ALL());
      auto Wi = Kokkos::subview(work, i, Kokkos::ALL(), Kokkos::ALL());
      Invert<3>(Fi, invFi);
      VoigtToMatrix(StressVoigti, Stressi);
      // W = Sigma@F^-T
      KokkosBatched::SerialGemm<
          KokkosBatched::Trans::NoTranspose,
          KokkosBatched::Trans::Transpose,
          KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, Stressi, invFi,
                                                        0.0, Wi);
      // Sigma (cauchy) = J F^-1 W
      auto J = Determinant<3>(Fi);
      KokkosBatched::SerialGemm<
          KokkosBatched::Trans::NoTranspose,
          KokkosBatched::Trans::NoTranspose,
          KokkosBatched::Algo::Gemm::Unblocked>::invoke(J, invFi, Wi, 0.0,
                                                        Stressi);
      // convert back to voigt
      MatrixToVoigt(Stressi, StressVoigti);
    });
  }
  // dF is increment in F
  template <typename T, typename... Args>
  auto compute_updated_deformation_gradient(
      Kokkos::View<T * [3][3], Args...> F,
      Kokkos::View<T * [3][3], Args...> dF) -> Kokkos::View<T * [3][3], Args...>
  {
    // TODO no initialize
    Kokkos::View<T * [3][3], Args...> F_updated("F", F.extent(0));
    Kokkos::parallel_for(
        F.extent(0), KOKKOS_LAMBDA(const int i) {
          using namespace KokkosBatched;
          using namespace Kokkos;
          auto dF_i = Kokkos::subview(dF, (dF.extent(0) > 1 ? i : 0),
                                      Kokkos::ALL(), Kokkos::ALL());
          auto F_i = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
          auto F_updated_i =
              Kokkos::subview(F_updated, i, Kokkos::ALL(), Kokkos::ALL());
          SerialGemm<Trans::NoTranspose, Trans::NoTranspose,
                     Algo::Gemm::Unblocked>::invoke(1.0, dF_i, F_i, 0.0,
                                                    F_updated_i);
        });
    return F_updated;
  }
  template <typename T, typename... Args>
  auto compute_left_cauchy_green(Kokkos::View<T * [3][3], Args...> F)
      -> Kokkos::View<T * [3][3], Args...>
  {
    // TODO no initialize
    Kokkos::View<T * [3][3], Args...> left_cauchy_green(
        "left cauchy green deformation tensor", F.extent(0));
    Kokkos::parallel_for(
        F.extent(0), KOKKOS_LAMBDA(const int i) {
          using namespace KokkosBatched;
          using namespace Kokkos;
          auto F_i = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
          auto left_cauchy_green_i = Kokkos::subview(
              left_cauchy_green, i, Kokkos::ALL(), Kokkos::ALL());
          SerialGemm<Trans::NoTranspose, Trans::Transpose,
                     Algo::Gemm::Unblocked>::invoke(1.0, F_i, F_i, 0.0,
                                                    left_cauchy_green_i);
        });
    return left_cauchy_green;
  }
  template<typename ViewT>
  auto FillIdentity(ViewT view) -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 3, void>
  {
    KOKKOS_ASSERT(view.extent(1) == 3);
    KOKKOS_ASSERT(view.extent(2) == 3);
      Kokkos::parallel_for(Kokkos::RangePolicy<typename ViewT::execution_space>(0, view.extent(0)),
        KOKKOS_LAMBDA(int i) {
          auto view_i = Kokkos::subview(view, i, Kokkos::ALL(), Kokkos::ALL());
          KokkosBatched::SerialSetIdentity::invoke(view_i);
        });

  }
  template <typename ViewT>
  KOKKOS_INLINE_FUNCTION auto ComputeGreenLagrangeStrain(ViewT F,
                                                         ViewT E) noexcept
  {
    KokkosBatched::SerialSetIdentity::invoke(E);
    KokkosBatched::SerialGemm<
        KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose,
        KokkosBatched::Algo::Gemm::Unblocked>::invoke(0.5, F, F, -0.5, E);
  }

}

#endif  // MUMFIM_SRC_MUMFIM_MICROSCALE_STRESSCONVERSION_H
