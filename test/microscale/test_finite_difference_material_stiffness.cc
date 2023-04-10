#include <mumfim/microscale/MicroTypeDefinitions.h>
#include <mumfim/microscale/MaterialStiffness.h>

#include <Catch2/catch.hpp>
#include <mumfim/microscale/BatchedNeohookeanAnalysis.h>
#include <mumfim/microscale/BatchedRVEAnalysis.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <chrono>

template <typename ViewT>
KOKKOS_INLINE_FUNCTION auto ComputeE(ViewT F, ViewT E) noexcept
{
  E(0, 0) = 1.0;
  E(1, 1) = 1.0;
  E(2, 2) = 1.0;
  KokkosBatched::SerialGemm<
      KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose,
      KokkosBatched::Algo::Gemm::Unblocked>::invoke(0.5, F, F, -0.5, E);
}

template <typename ExeSpace>
struct LinearStress
{
  using memory_space = typename ExeSpace::memory_space;
  auto operator()(
      Kokkos::View<double * [3][3], memory_space> F,
      Kokkos::View<double * [3][3], memory_space> /*unused*/) noexcept
      -> Kokkos::View<double * [6], memory_space>
  {
    Kokkos::View<double * [3][3], typename ExeSpace::memory_space> E_matrix(
        "E_matrix", F.extent(0));
    Kokkos::View<double * [6], typename ExeSpace::memory_space> E("E",
                                                                  F.extent(0));
    Kokkos::deep_copy(E_matrix, 0.0);
    Kokkos::deep_copy(E, 0.0);
    Kokkos::parallel_for(
        "linear stress", Kokkos::RangePolicy<ExeSpace>(0, F.extent(0)),
        KOKKOS_LAMBDA(int i) {
          auto F_i = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
          auto E_matrix_i =
              Kokkos::subview(E_matrix, i, Kokkos::ALL(), Kokkos::ALL());
          auto E_i = Kokkos::subview(E, i, Kokkos::ALL());
          ComputeE(F_i, E_matrix_i);
          mumfim::MatrixToVoigt(E_matrix_i, E_i);
        });
    Kokkos::fence();
    return E;
  }
};

template <typename ExeSpace>
struct CubicStress
{
  using memory_space = typename ExeSpace::memory_space;
  auto operator()(Kokkos::View<double * [3][3], memory_space> F,
                  Kokkos::View<double * [3][3], memory_space> /*unused*/,
                  bool /*unused*/ = false) noexcept
      -> Kokkos::View<double * [6], memory_space>
  {
    Kokkos::View<double * [3][3], typename ExeSpace::memory_space> E_matrix(
        "E_matrix", F.extent(0));
    Kokkos::View<double * [3][3], typename ExeSpace::memory_space> tmp(
        "tmp", F.extent(0));
    Kokkos::View<double * [3][3], typename ExeSpace::memory_space> tmp2(
        "tmp", F.extent(0));
    Kokkos::View<double * [6], typename ExeSpace::memory_space> E("E",
                                                                  F.extent(0));
    Kokkos::deep_copy(E_matrix, 0.0);
    Kokkos::deep_copy(E, 0.0);
    Kokkos::parallel_for(
        "linear stress", Kokkos::RangePolicy<ExeSpace>(0, F.extent(0)),
        KOKKOS_LAMBDA(int i) {
          auto F_i = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
          auto E_matrix_i =
              Kokkos::subview(E_matrix, i, Kokkos::ALL(), Kokkos::ALL());
          auto tmpi = Kokkos::subview(tmp, i, Kokkos::ALL(), Kokkos::ALL());
          auto tmp2i = Kokkos::subview(tmp2, i, Kokkos::ALL(), Kokkos::ALL());
          auto E_i = Kokkos::subview(E, i, Kokkos::ALL());
          ComputeE(F_i, E_matrix_i);
          KokkosBatched::SerialGemm<
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, E_matrix_i,
                                                            E_matrix_i, 0.0,
                                                            tmpi);
          KokkosBatched::SerialGemm<
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, E_matrix_i,
                                                            tmpi, 0.0, tmp2i);
          mumfim::MatrixToVoigt(tmp2i, E_i);
        });
    Kokkos::fence();
    return E;
  }
};

TEST_CASE("dUdE", "[finite_difference]")
{
  Kokkos::Random_XorShift64_Pool<Kokkos::HostSpace> random(13718);
  constexpr int num_tests = 10;
  Kokkos::View<double * [3][3], Kokkos::HostSpace> F("F", num_tests);
  Kokkos::View<double * [6], Kokkos::HostSpace> PK2("PK2", num_tests);
  Kokkos::fill_random(F, random, 1.0);
  Kokkos::parallel_for(
      F.extent(0), KOKKOS_LAMBDA(int i) {
        auto Fi = Kokkos::subview(F,i, Kokkos::ALL(), Kokkos::ALL());
        while(mumfim::Determinant<3>(Fi) < 0.0)
        {
          Fi(0, 0) += 1.0;
          Fi(1, 1) += 1.0;
          Fi(2, 2) += 1.0;
        }
      });
  Kokkos::parallel_for(
      F.extent(0), KOKKOS_LAMBDA(int i) {
        auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
        KOKKOS_ASSERT(mumfim::Determinant<3>(Fi) > 0.0);
      });

  auto t1 = std::chrono::steady_clock::now();
  // mumfim::computeDPK2dE<Kokkos::Serial>(F, dPK2dU);
  mumfim::BatchedNeohookeanAnalysis analysis(F.extent(0), 1E3, 0.3);
  using exe_space = typename decltype(analysis)::exe_space;
  mumfim::BatchedAnalysisGetPK2StressFunc compute_pk2_stress{analysis};
  //  LinearStress<exe_space> compute_pk2_stress{};
  // CubicStress<exe_space> compute_pk2_stress{};
  auto [U, R] = mumfim::PolarDecomposition<exe_space>(F);
  // apply the full deformation gradient as the increment
  // in the first timestep
  auto stress = compute_pk2_stress(F, F, true);
  //mumfim::StressFiniteDifferenceFunc finite_difference{compute_pk2_stress,
  //                                                     stress, 1E-7};
  mumfim::StressCentralDifferenceFunc finite_difference{compute_pk2_stress, 1E-7};
  auto result = mumfim::ComputeDPK2dE<exe_space>(F, finite_difference);
  mumfim::ConvertTLStiffnessToULStiffness<exe_space>(F, result);
  Kokkos::DualView<double * [6][6], typename exe_space::memory_space> stiffness(
      "stiffness", F.extent(0));
  analysis.computeMaterialStiffness(stiffness);
  stiffness.sync_host();
  auto stiffness_h = stiffness.h_view;
  for (int l=0; l<num_tests; ++l) {

  for (size_t i = 0; i < result.extent(1); ++i)
  {
    for (size_t j = 0; j < result.extent(2); ++j)
    {
      auto r = result(l, i, j) < 1E-5 ? 0.0 : result(l, i, j);
      std::cout << r << " ";
    }
    std::cout << "\n";
  }
  std::cout << "---------------------\n";
  for (size_t i = 0; i < stiffness_h.extent(1); ++i)
  {
    for (size_t j = 0; j < stiffness_h.extent(2); ++j)
    {
      auto r = stiffness_h(l, i, j) < 1E-5 ? 0.0 : stiffness_h(l, i, j);
      std::cout << r << " ";
    }
    std::cout << "\n";
  }
  }
  auto t2 = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration<double>(t2 - t1).count();
  std::cout << "Time taken by function: " << duration << " seconds"
            << std::endl;
}
