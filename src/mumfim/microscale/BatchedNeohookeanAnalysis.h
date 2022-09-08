#ifndef MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDNEOHOOKEANANALYSIS_H
#define MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDNEOHOOKEANANALYSIS_H
#include <Kokkos_core.hpp>
#include "BatchedRVEAnalysis.h"
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "MicroTypeDefinitions.h"
namespace mumfim
{
  template <typename T, typename... Args>
  KOKKOS_INLINE_FUNCTION auto determinant(Kokkos::View<T[3][3], Args...> matrix)
  {
    auto a = matrix(0, 0), b = matrix(0, 1), c = matrix(0, 2);
    auto d = matrix(1, 0), e = matrix(1, 1), f = matrix(1, 2);
    auto g = matrix(1, 0), h = matrix(1, 1), i = matrix(1, 2);
    return a * e * i - a * f * h + b * f * g - b * d * i + c * d * h -
           c * e * g;
  }
  template <typename Scalar = mumfim::Scalar,
            typename LocalOrdinal = mumfim::LocalOrdinal,
            typename ExeSpace = Kokkos::DefaultExecutionSpace>
  class BatchedNeohookeanAnalysis
      : public BatchedRVEAnalysis<Scalar, LocalOrdinal, ExeSpace>
  {
    public:
    explicit BatchedNeohookeanAnalysis(LocalOrdinal num_rves)
        : BatchedRVEAnalysis<Scalar, LocalOrdinal, ExeSpace>(num_rves)
        , F("f", num_rves)
    {
      // initialize original deformation gradient
      Kokkos::parallel_for(
          F.extent(0), KOKKOS_LAMBDA(int i) {
            for (int j = 0; j < 3; ++j)
            {
              F(i, j, j) = 1;
            }
          });
    }
    bool run(Kokkos::DualView<Scalar * [3][3], ExeSpace> dfmGrds,
             Kokkos::DualView<Scalar * [6], ExeSpace> sigma,
             bool update_coords = true) final
    {
      dfmGrds.sync_device();
      dfmGrds.sync_host();
      printf("dfmgrds: %f %f %f %f %f %f %f %f %f\n", dfmGrds.h_view(0, 0, 0),
             dfmGrds.h_view(0, 0, 1), dfmGrds.h_view(0, 0, 2),
             dfmGrds.h_view(0, 1, 0), dfmGrds.h_view(0, 1, 1),
             dfmGrds.h_view(0, 1, 2), dfmGrds.h_view(0, 2, 0),
             dfmGrds.h_view(0, 2, 1), dfmGrds.h_view(0, 2, 2));
      Kokkos::View<Scalar * [3][3]> green_strain("green_strain", F.extent(0));
      Kokkos::View<Scalar * [3][3]> F_updated("F_updated", F.extent(0));
      Kokkos::parallel_for(
          F.extent(0), KOKKOS_LAMBDA(const int i) {
            using namespace KokkosBatched;
            using namespace Kokkos;
            auto Fi_incremental =
                Kokkos::subview(dfmGrds.d_view, (dfmGrds.extent(0) > 1 ? i : 0),
                                Kokkos::ALL(), Kokkos::ALL());
            auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
            auto F_updated_i =
                Kokkos::subview(F_updated, i, Kokkos::ALL(), Kokkos::ALL());
            if (i == 0)
              printf("I%f, %f, %f\n", Fi_incremental(0, 0),
                     Fi_incremental(1, 1), Fi_incremental(2, 2));
            if (i == 0) printf("A%f, %f, %f\n", Fi(0, 0), Fi(1, 1), Fi(2, 2));
            auto sigma_i = Kokkos::subview(sigma.d_view, i, Kokkos::ALL());
            auto green_strain_i =
                Kokkos::subview(green_strain, i, Kokkos::ALL(), Kokkos::ALL());
            SerialGemm<Trans::NoTranspose, Trans::NoTranspose,
                       Algo::Gemm::Unblocked>::invoke(1.0, Fi_incremental, Fi,
                                                      0.0, F_updated_i);
            // printf("B%f, %f, %f\n", Fi(0,0), Fi(1,1), Fi(2,2));
            if (i == 0)
              printf("B%f, %f, %f\n", F_updated_i(0, 0), F_updated_i(1, 1),
                     F_updated_i(2, 2));
            SerialGemm<Trans::Transpose, Trans::NoTranspose,
                       Algo::Gemm::Unblocked>::invoke(1.0 / 2.0, F_updated_i,
                                                      F_updated_i, 0.0,
                                                      green_strain_i);
            green_strain_i(0, 0) -= 0.5;
            green_strain_i(1, 1) -= 0.5;
            green_strain_i(2, 2) -= 0.5;
            double prefact = youngs_modulus /
                             ((1 + poissons_ratio) * (1 - 2 * poissons_ratio));
            sigma_i(0) =
                prefact * (green_strain_i(0, 0) * (1 - poissons_ratio) +
                           green_strain_i(1, 1) * poissons_ratio +
                           green_strain_i(2, 2) * poissons_ratio);
            sigma_i(1) =
                prefact * (green_strain_i(1, 1) * (1 - poissons_ratio) +
                           green_strain_i(0, 0) * poissons_ratio +
                           green_strain_i(2, 2) * poissons_ratio);
            sigma_i(2) =
                prefact * (green_strain_i(2, 2) * (1 - poissons_ratio) +
                           green_strain_i(1, 1) * poissons_ratio +
                           green_strain_i(0, 0) * poissons_ratio);
            sigma_i(3) = prefact * (green_strain_i(1, 2) *
                                    (1.0 - 2 * poissons_ratio) / 2);
            sigma_i(4) = prefact * (green_strain_i(0, 2) *
                                    (1.0 - 2 * poissons_ratio) / 2);
            sigma_i(5) = prefact * (green_strain_i(0, 1) *
                                    (1.0 - 2 * poissons_ratio) / 2);
            if (i == 0)
              printf("Sigma: %f, %f, %f, %f, %f, %f\n", sigma_i(0), sigma_i(1),
                     sigma_i(2), sigma_i(3), sigma_i(4), sigma_i(5));
            if (i == 0)
              printf("GS: %f:%f, %f:%f, %f:%f, %f:%f, %f:%f, %f:%f\n",
                     green_strain_i(0, 0), pow(Fi_incremental(0, 0), 2),
                     green_strain_i(1, 1), pow(Fi_incremental(1, 1), 2),
                     green_strain_i(2, 2), pow(Fi_incremental(2, 2), 2),
                     green_strain_i(1, 2), pow(Fi_incremental(1, 2), 2),
                     green_strain_i(0, 2), pow(Fi_incremental(0, 2), 2),
                     green_strain_i(0, 1), pow(Fi_incremental(0, 1), 2));
          });
      sigma.modify_device();
      if (update_coords)
      {
        Kokkos::deep_copy(F, F_updated);
      }
      return true;
    }
    void computeMaterialStiffness(
        Kokkos::DualView<Scalar * [6][6], ExeSpace> C) final
    {
      C.sync_device();
      Kokkos::parallel_for(
          C.extent(0), KOKKOS_LAMBDA(int i) {
            auto C_i =
                Kokkos::subview(C.d_view, i, Kokkos::ALL(), Kokkos::ALL());
            C_i(0, 0) = (1 / youngs_modulus) * 1;
            C_i(1, 1) = (1 / youngs_modulus) * 1;
            C_i(2, 2) = (1 / youngs_modulus) * 1;
            C_i(3, 3) = (1 / youngs_modulus) * 2 * (1 + poissons_ratio);
            C_i(4, 4) = (1 / youngs_modulus) * 2 * (1 + poissons_ratio);
            C_i(5, 5) = (1 / youngs_modulus) * 2 * (1 + poissons_ratio);
            C_i(0, 1) = C_i(1, 0) = (1 / youngs_modulus) * -poissons_ratio;
            C_i(0, 2) = C_i(2, 0) = (1 / youngs_modulus) * -poissons_ratio;
          });
    }

    private:
    Kokkos::View<Scalar * [3][3], ExeSpace> F;
    static constexpr double youngs_modulus{0.0001};
    static constexpr double poissons_ratio{0.3};
  };
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDNEOHOOKEANANALYSIS_H
