#ifndef MUMFIM_SRC_MUMFIM_MICROSCALE_STRESSCONVERSION_H
#define MUMFIM_SRC_MUMFIM_MICROSCALE_STRESSCONVERSION_H
#include <Kokkos_Core.hpp>
#include "mumfim/microscale/MicroTypeDefinitions.h"
#include "mumfim/microscale/TensorUtilities.h"
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>

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

}

#endif  // MUMFIM_SRC_MUMFIM_MICROSCALE_STRESSCONVERSION_H
