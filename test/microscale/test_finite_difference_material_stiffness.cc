#include <mumfim/microscale/BatchedNeohookeanAnalysis.h>
#include <mumfim/microscale/BatchedRVEAnalysis.h>
#include <mumfim/microscale/MaterialStiffness.h>
#include <mumfim/microscale/MicroTypeDefinitions.h>

#include <catch2/catch.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

template <typename ViewT>
KOKKOS_INLINE_FUNCTION auto ComputeGreenLagrangeStrain(ViewT F,
                                                       ViewT E) noexcept
{
  E(0, 0) = 1.0;
  E(1, 1) = 1.0;
  E(2, 2) = 1.0;
  KokkosBatched::SerialGemm<
      KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose,
      KokkosBatched::Algo::Gemm::Unblocked>::invoke(0.5, F, F, -0.5, E);
}

template <typename ExeSpace>
struct LinearU
{
  using memory_space = typename ExeSpace::memory_space;
  using exe_space = ExeSpace;
  auto operator()(
      Kokkos::View<double * [3][3], memory_space> F,
      Kokkos::View<double * [3][3], memory_space> /*unused*/) noexcept
      -> Kokkos::View<double * [6], memory_space>
  {
    auto polar_result = mumfim::PolarDecomposition<exe_space>(F);
    KOKKOS_ASSERT(F.extent(0) == polar_result.U.extent(0));
    KOKKOS_ASSERT(F.extent(0) == polar_result.R.extent(0));
    Kokkos::View<double * [6], typename ExeSpace::memory_space> U_voigt(
        "U", F.extent(0));
    Kokkos::parallel_for(
        "linear stress", Kokkos::RangePolicy<ExeSpace>(0, F.extent(0)),
        KOKKOS_LAMBDA(int i) {
          auto U_i = Kokkos::subview(polar_result.U, i, Kokkos::ALL(), Kokkos::ALL());
          auto U_voigt_i = Kokkos::subview(U_voigt, i, Kokkos::ALL());
          mumfim::MatrixToVoigt(U_i, U_voigt_i);
        });
    return U_voigt;
  }
};

template <typename ExeSpace>
struct LinearStress
{
  using memory_space = typename ExeSpace::memory_space;
  using exe_space = ExeSpace;
  auto operator()(Kokkos::View<double * [3][3], memory_space> F,
                  Kokkos::View<double * [3][3], memory_space> /*unused*/,
                  bool /*unused*/ = false) noexcept
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
          ComputeGreenLagrangeStrain(F_i, E_matrix_i);
          mumfim::MatrixToVoigt(E_matrix_i, E_i);
        });
    return E;
  }
};

template <typename ExeSpace>
struct CubicStress
{
  using exe_space = ExeSpace;
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
          ComputeGreenLagrangeStrain(F_i, E_matrix_i);
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
    return E;
  }
};

template <typename ViewT, typename ViewT2>
[[nodiscard]] auto ViewsEqual(ViewT a,
                              ViewT2 b,
                              double rtol = 1E-9,
                              double atol = 0.0)
    -> std::enable_if_t<Kokkos::is_view<ViewT>::value &&
                            Kokkos::is_view<ViewT2>::value &&
                            ViewT::rank == 3 && ViewT2::rank == 3,
                        int>
{
  if (a.extent(0) != b.extent(0))
  {
    return false;
  }
  if (a.extent(1) != b.extent(1))
  {
    return false;
  }
  if (a.extent(2) != b.extent(2))
  {
    return false;
  }
  int sum = 0;
  Kokkos::parallel_reduce(
      "views equal",
      Kokkos::RangePolicy<typename ViewT::execution_space>(0, a.extent(0)),
      KOKKOS_LAMBDA(int i, int & update) {
        for (size_t j = 0; j < a.extent(1); ++j)
        {
          for (size_t k = 0; k < a.extent(2); ++k)
          {
            bool is_close =
                (fabs(a(i, j, k) - b(i, j, k)) <=
                 Kokkos::max(rtol * Kokkos::max(fabs(a(i, j, k)), fabs(b(i, j, k))),
                          atol));
            update += !is_close;
          }
        }
      },
      sum);
  return (sum == 0);
}

// MandelToVoight -> VoigtToMandel gives the same matrix

TEST_CASE("finite_difference", "[material_stiffness],[finite_difference]")
{
  using exe_space = Kokkos::DefaultExecutionSpace;
  //using exe_space = Kokkos::Serial;
  using memory_space = typename exe_space::memory_space;
  // finite difference should return dPK2dU. So, we can analyze two cases.
  Kokkos::Random_XorShift64_Pool<memory_space> random(13718);
  constexpr int num_tests = 10;
  Kokkos::View<double * [3][3], memory_space> F("F", num_tests);
  Kokkos::View<double * [6], memory_space> PK2("PK2", num_tests);
  Kokkos::fill_random(F, random, 1.0);
  Kokkos::parallel_for(Kokkos::RangePolicy<exe_space>(0, F.extent(0)), KOKKOS_LAMBDA(int i) {
        auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
        while (mumfim::Determinant<3>(Fi) < 0.5)
        {
          Fi(0, 0) += 1.0;
          Fi(1, 1) += 1.0;
          Fi(2, 2) += 1.0;
        }
      });
  Kokkos::parallel_for(Kokkos::RangePolicy<exe_space>(0, F.extent(0)), KOKKOS_LAMBDA(int i) {
        auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
        KOKKOS_ASSERT(mumfim::Determinant<3>(Fi) > 0.0);
      });
  Kokkos::View<double * [6], memory_space> current_stress("current_stress",
                                                          num_tests);
  auto [U, R] = mumfim::PolarDecomposition<exe_space>(F);
  auto F_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, F);
  auto R_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, R);
  auto U_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, U);

  auto T = mumfim::ComputeT<memory_space>();
  auto T_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, T);
  // 1. Stress function is linear in U i.e., the resulting derivative should be
  // the identity matrix
  SECTION("Linear Stress")
  {
    LinearU<exe_space> compute_pk2_stress{};
    current_stress = compute_pk2_stress(F, F);
    auto current_stress_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, current_stress);
    SECTION("Forward Difference")
    {
      mumfim::StressFiniteDifferenceFunc compute_dpk2du(compute_pk2_stress,
                                                        current_stress);
      auto dPK2dU = compute_dpk2du(F, R, U);
      auto dPK2dU_h =
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dPK2dU);
      for (int i = 0; i < num_tests; ++i)
      {
        for (int j = 0; j < 6; ++j)
        {
          for (int k = 0; k < 6; ++k)
          {
            REQUIRE(dPK2dU_h(i, j, k) == Approx(T_h(j, k)).margin(1E-6));
          }
        }
      }
    }
    SECTION("Central Difference")
    {
      mumfim::StressCentralDifferenceFunc compute_dpk2du(compute_pk2_stress);
      auto dPK2dU = compute_dpk2du(F, R, U);
      auto dPK2dU_h =
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dPK2dU);
      for (int i = 0; i < num_tests; ++i)
      {
        for (int j = 0; j < 6; ++j)
        {
          for (int k = 0; k < 6; ++k)
          {
            REQUIRE(dPK2dU_h(i, j, k) == Approx(T_h(j, k)).margin(1E-6));
          }
        }
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION
int Delta(int i, int j) { return (i == j) ? 1 : 0; }

// 1. Get PK2 stress
// 2. Finite Difference
// 3. dPK2/dE
TEST_CASE("dPK2dE", "[material stiffness]")
{
  using exe_space = Kokkos::DefaultExecutionSpace;
  using memory_space = typename Kokkos::DefaultExecutionSpace::memory_space;
  Kokkos::Random_XorShift64_Pool<memory_space> random(13718);
  constexpr int num_tests = 100;
  Kokkos::View<double * [3][3], memory_space> F("F", num_tests);
  Kokkos::View<double * [6], memory_space> PK2("PK2", num_tests);
  Kokkos::fill_random(F, random, 1.0);
  Kokkos::parallel_for(
      F.extent(0), KOKKOS_LAMBDA(int i) {
        auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
        while (mumfim::Determinant<3>(Fi) < 0.5)
        {
          Fi(0, 0) += 0.2;
          Fi(1, 1) += 0.1;
          Fi(2, 2) += 0.1;
        }
      });
  Kokkos::parallel_for(
      F.extent(0), KOKKOS_LAMBDA(int i) {
        auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
        KOKKOS_ASSERT(mumfim::Determinant<3>(Fi) > 0.0);
      });
  SECTION("Linear Stress")
  {
    LinearStress<exe_space> compute_pk2_stress{};
    // apply the full deformation gradient as the increment
    // in the first timestep
    auto stress = compute_pk2_stress(F, F, true);
    Kokkos::View<double * [6][6], memory_space> dPK2dU_reference(
        "dPK2dU_reference", num_tests);
    Kokkos::parallel_for(
        "set reference dpk2dU", Kokkos::RangePolicy<exe_space>(0, num_tests),
        KOKKOS_LAMBDA(int i) {
          for (int j = 0; j < 6; ++j)
          {
            for (int k = 0; k < 6; ++k)
            {
              if (j == k)
              {
                if (j < 3)
                {
                  dPK2dU_reference(i, j, k) = 1.0;
                }
                else
                {
                  dPK2dU_reference(i, j, k) = 0.5;
                }
              }
            }
          }
        });
    SECTION("Forward Difference")
    {
      mumfim::StressFiniteDifferenceFunc dPK2dU_func(compute_pk2_stress, stress,
                                                     1E-7);
      auto dPK2dE = mumfim::ComputeDPK2dE<exe_space>(F, dPK2dU_func);
      REQUIRE(ViewsEqual(dPK2dE, dPK2dU_reference, 1E-8, 1E-4));
    }
    SECTION("Central Difference")
    {
      mumfim::StressCentralDifferenceFunc dPK2dU_func(compute_pk2_stress, 1E-7);
      auto dPK2dE = mumfim::ComputeDPK2dE<exe_space>(F, dPK2dU_func);
      REQUIRE(ViewsEqual(dPK2dE, dPK2dU_reference, 1E-8, 1E-8));
    }
  }
  SECTION("Cubic Stress")
  {
    CubicStress<exe_space> compute_pk2_stress{};
    // apply the full deformation gradient as the increment
    // in the first timestep
    auto stress = compute_pk2_stress(F, F, true);
    Kokkos::View<double * [6][6], memory_space> dPK2dU_reference(
        "dPK2dU_reference", num_tests);
    Kokkos::deep_copy(dPK2dU_reference, 0.0);
    Kokkos::View<double * [3][3], memory_space> E("E", num_tests);

    Kokkos::parallel_for(
        "compute GL strain", Kokkos::RangePolicy<exe_space>(0, num_tests),
        KOKKOS_LAMBDA(int i) {
          auto Ei = Kokkos::subview(E, i, Kokkos::ALL(), Kokkos::ALL());
          auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
          ComputeGreenLagrangeStrain(Fi, Ei);
        });
    Kokkos::parallel_for(
        "set reference dpk2dU", Kokkos::RangePolicy<exe_space>(0, num_tests),
        KOKKOS_LAMBDA(int p) {
          for (int i = 0; i < 3; ++i)
          {
            for (int j = i; j < 3; ++j)
            {
              for (int k = 0; k < 3; ++k)
              {
                for (int l = k; l < 3; ++l)
                {
                  double sum = 0.0;
                  // sum indicies: q
                  sum += E(p, i, k) * E(p, l, j);
                  sum += E(p, l, i) * E(p, j, k);
                  for (int q = 0; q < 3; ++q)
                  {
                    sum += Delta(i, k) * E(p, l, q) * E(p, q, j);
                    sum += Delta(j, l) * E(p, i, q) * E(p, q, k);
                    sum += Delta(i, l) * E(p, q, k) * E(p, j, q);
                    sum += Delta(j, k) * E(p, q, i) * E(p, l, q);
                  }
                  dPK2dU_reference(p, mumfim::TensorIndex2VoigtIndex(i, j),
                                   mumfim::TensorIndex2VoigtIndex(k, l)) +=
                      (sum / 2.0);
                }
              }
            }
          }
        });
    SECTION("Forward Difference")
    {
      mumfim::StressFiniteDifferenceFunc dPK2dU_func(compute_pk2_stress, stress,
                                                     1E-7);
      auto dPK2dE = mumfim::ComputeDPK2dE<exe_space>(F, dPK2dU_func);
      REQUIRE(ViewsEqual(dPK2dE, dPK2dU_reference, 1E-5, 1E-5));
    }
    SECTION("Central Difference")
    {
      mumfim::StressCentralDifferenceFunc dPK2dU_func(compute_pk2_stress, 1E-7);
      auto dPK2dE = mumfim::ComputeDPK2dE<exe_space>(F, dPK2dU_func);
      REQUIRE(ViewsEqual(dPK2dE, dPK2dU_reference, 1E-6, 1E-6));
    }
  }
  SECTION("NeoHookean")
  {
    SECTION("Forward Difference")
    {
      mumfim::BatchedNeohookeanAnalysis analysis(F.extent(0), 1E3, 0.3);
      static_assert(
          std::is_same_v<typename decltype(analysis)::exe_space, exe_space>);
      mumfim::BatchedAnalysisGetPK2StressFunc compute_pk2_stress{analysis};
      // apply the full deformation gradient as the increment
      // in the first timestep
      auto stress = compute_pk2_stress(F, F, true);
      mumfim::StressFiniteDifferenceFunc dPK2dU_func(compute_pk2_stress, stress,
                                                     1E-7);
      auto dPK2dE = mumfim::ComputeDPK2dE<exe_space>(F, dPK2dU_func);
      mumfim::ConvertTLStiffnessToULStiffness<exe_space>(F, dPK2dE);

      // the analytic stress function
      Kokkos::DualView<double * [6][6], typename exe_space::memory_space>
          stiffness("stiffness", F.extent(0));
      analysis.computeMaterialStiffness(stiffness);
      stiffness.sync_device();
      REQUIRE(ViewsEqual(stiffness.d_view, dPK2dE, 1E-8, 1E-2));
    }
    SECTION("Central Difference")
    {
      mumfim::BatchedNeohookeanAnalysis analysis(F.extent(0), 1E3, 0.3);
      static_assert(
          std::is_same_v<typename decltype(analysis)::exe_space, exe_space>);
      mumfim::BatchedAnalysisGetPK2StressFunc compute_pk2_stress{analysis};
      // apply the full deformation gradient as the increment
      // in the first timestep
      compute_pk2_stress(F, F, true);
      mumfim::StressCentralDifferenceFunc dPK2dU_func(compute_pk2_stress, 1E-7);
      auto dPK2dE = mumfim::ComputeDPK2dE<exe_space>(F, dPK2dU_func);
      mumfim::ConvertTLStiffnessToULStiffness<exe_space>(F, dPK2dE);


      // the analytic stress function
      Kokkos::DualView<double * [6][6], typename exe_space::memory_space>
          stiffness("stiffness", F.extent(0));
      analysis.computeMaterialStiffness(stiffness);
      stiffness.sync_device();
      stiffness.sync_host();

      auto dPK2dE_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dPK2dE);
      for(int i=0; i<dPK2dE_h.extent_int(0); ++i) {
        for(int j=0;j<dPK2dE_h.extent_int(1); ++j) {
          for(int k=0;k<dPK2dE_h.extent_int(2); ++k) {
            std::cout << dPK2dE_h(i,j,k) << " ";
          }
          std::cout<<"\n";
        }
        std::cout<<"----------------\n";
        for(int j=0;j<stiffness.h_view.extent_int(1); ++j) {
          for(int k=0;k<stiffness.h_view.extent_int(2); ++k) {
            std::cout << stiffness.h_view(i,j,k) << " ";
          }
          std::cout<<"\n";
        }
        std::cout<<"=====================\n";
      }
      REQUIRE(ViewsEqual(stiffness.d_view, dPK2dE, 1E-8, 1E-3));
    }
  }
}
