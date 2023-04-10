#ifndef MUMFIM_SRC_MUMFIM_MICROSCALE_MATERIALSTIFFNESS_H
#define MUMFIM_SRC_MUMFIM_MICROSCALE_MATERIALSTIFFNESS_H
#include <mumfim/microscale/MicroTypeDefinitions.h>
#include <mumfim/microscale/StressConversion.h>
#include <mumfim/microscale/TensorUtilities.h>
#include <mumfim/microscale/PolarDecomposition.h>
#include <KokkosBatched_Scale_Decl.hpp>
#include <KokkosBatched_Scale_Impl.hpp>
#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Copy_Impl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
#include <KokkosBatched_SolveLU_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Team_Impl.hpp>

#include <Kokkos_Core.hpp>

namespace mumfim
{
  template <typename ExeSpace>
  auto ComputeM(
      Kokkos::View<double * [3][3], typename ExeSpace::memory_space> U)
      -> Kokkos::View<double * [6][6], typename ExeSpace::memory_space>
  {
    Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space> M(
        "M", U.extent(0));
    Kokkos::deep_copy(M, 0.0);
    Kokkos::parallel_for(
        "FillM", Kokkos::RangePolicy<ExeSpace>(0, U.extent(0)),
        KOKKOS_LAMBDA(int i) {
          auto sr2 = sqrt(2.0);
          auto Mi = Kokkos::subview(M, i, Kokkos::ALL(), Kokkos::ALL());
          Mi(0, 0) = 2 * U(i, 0, 0);
          Mi(0, 4) = sr2 * U(i, 0, 2);
          Mi(0, 5) = sr2 * U(i, 0, 1);
          Mi(1, 1) = 2 * U(i, 1, 1);
          Mi(1, 3) = sr2 * U(i, 1, 2);
          Mi(1, 5) = sr2 * U(i, 0, 1);
          Mi(2, 2) = 2 * U(i, 2, 2);
          Mi(2, 3) = sr2 * U(i, 1, 2);
          Mi(2, 4) = sr2 * U(i, 0, 2);
          Mi(3, 1) = sr2 * U(i, 1, 2);
          Mi(3, 2) = sr2 * U(i, 1, 2);
          Mi(3, 3) = U(i, 1, 1) + U(i, 2, 2);
          Mi(3, 4) = U(i, 0, 1);
          Mi(3, 5) = U(i, 0, 2);
          Mi(4, 0) = sr2 * U(i, 0, 2);
          Mi(4, 2) = sr2 * U(i, 0, 2);
          Mi(4, 3) = U(i, 0, 1);
          Mi(4, 4) = U(i, 0, 0) + U(i, 2, 2);
          Mi(4, 5) = U(i, 1, 2);
          Mi(5, 0) = sr2 * U(i, 0, 1);
          Mi(5, 1) = sr2 * U(i, 0, 1);
          Mi(5, 3) = U(i, 0, 2);
          Mi(5, 4) = U(i, 1, 2);
          Mi(5, 5) = U(i, 0, 0) + U(i, 1, 1);
          KokkosBatched::SerialScale::invoke(0.5, Mi);
        });
    return M;
  }
  // calculate the mandel form of the probing directions
  template <typename MemorySpace>
  auto ComputeT() -> Kokkos::View<Scalar[6][6], MemorySpace>
  {
    Kokkos::View<Scalar[6][6], MemorySpace> T("T");
    auto T_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, T);
    Kokkos::deep_copy(T_h, 0.0);
    auto sr2 = sqrt(2);
    T_h(0, 0) = 1;
    T_h(1, 1) = 1;
    T_h(2, 2) = 1;
    T_h(3, 3) = 2;
    T_h(4, 4) = 2;
    T_h(5, 5) = 2;
    T_h(0, 4) = sr2;
    T_h(0, 5) = sr2;
    T_h(1, 3) = sr2;
    T_h(1, 5) = sr2;
    T_h(2, 3) = sr2;
    T_h(2, 4) = sr2;
    Kokkos::deep_copy(T, T_h);
    return T;
  }

  // These probe vectors correspond to the mandel form that is created in
  // the ComputeT function
  template <typename MemorySpace>
  auto ComputeProbeVectors() -> Kokkos::View<Scalar[6][3][3], MemorySpace>
  {
    Kokkos::View<Scalar[6][3][3], MemorySpace> T("T");
    auto T_h = Kokkos::create_mirror_view(Kokkos::HostSpace{}, T);
    Kokkos::deep_copy(T_h, 0.0);
    T_h(0, 0, 0) = 1.0;
    T_h(1, 1, 1) = 1.0;
    T_h(2, 2, 2) = 1.0;
    T_h(3, 1, 1) = 1.0;
    T_h(3, 1, 2) = 1.0;
    T_h(3, 2, 1) = 1.0;
    T_h(3, 2, 2) = 1.0;
    T_h(4, 0, 0) = 1.0;
    T_h(4, 0, 2) = 1.0;
    T_h(4, 2, 0) = 1.0;
    T_h(4, 2, 2) = 1.0;
    T_h(5, 0, 0) = 1.0;
    T_h(5, 0, 1) = 1.0;
    T_h(5, 1, 0) = 1.0;
    T_h(5, 1, 1) = 1.0;
    Kokkos::deep_copy(T, T_h);
    return T;
  }

  /**
   * @brief Compute the derivative of the PK2 stress with respect to the Green
   * Lagrange Strain This is the material stiffness that is used in the
   * Total-Lagrangian formulation
   * @tparam ExeSpace Execution space to use
   * @tparam Func (Kokkos::View<Scalar*[3][3]> U) -> Kokkos::View<Scalar*[6]>
   * (dPk2_dU)
   * @param F deformation gradient
   * @param compute_dpk2_dU function to compute the derivative of the PK2 stress
   * with respect to the stretch tensor
   * @return derivative of the PK2 stress with respect to the Green Lagrange
   * Strain (Total Lagrangian Material Stiffness) in VoigtForm
   */
  template <typename ExeSpace, typename Func>
  auto ComputeDPK2dE(
      Kokkos::View<Scalar * [3][3], typename ExeSpace::memory_space> F,
      Func & compute_dpk2_dU)
      -> Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space>
  {
    using memory_space = typename ExeSpace::memory_space;
    auto [U, R] = PolarDecomposition<ExeSpace>(F);
    auto M = ComputeM<ExeSpace>(U);
    auto T = ComputeT<memory_space>();

    // auto P = ComputeDiff<ExeSpace>()
    //  TODO COMPUTE P
    Kokkos::View<Scalar * [6][6], memory_space> P = compute_dpk2_dU(F, R, U);
    Kokkos::fence();
    Kokkos::View<Scalar * [6][6], memory_space> MT("MT", U.extent(0));
    Kokkos::View<Scalar * [6][6], memory_space> D("D", U.extent(0));
    Kokkos::View<Scalar * [6][6], memory_space> PTr("P Transpose", U.extent(0));
    auto policy = Kokkos::TeamPolicy<ExeSpace>(U.extent(0), Kokkos::AUTO());
    using member_type = typename decltype(policy)::member_type;
    Kokkos::parallel_for(
        "compute MT", policy, KOKKOS_LAMBDA(member_type member) {
          auto i = member.league_rank();
          // D^T = Solve( ( (MT@T) ^T, P^T)
          auto Mi = Kokkos::subview(M, i, Kokkos::ALL(), Kokkos::ALL());
          auto MTi = Kokkos::subview(MT, i, Kokkos::ALL(), Kokkos::ALL());
          auto Pi = Kokkos::subview(P, i, Kokkos::ALL(), Kokkos::ALL());
          auto PTri = Kokkos::subview(PTr, i, Kokkos::ALL(), Kokkos::ALL());
          auto Di = Kokkos::subview(D, i, Kokkos::ALL(), Kokkos::ALL());
          KokkosBatched::TeamCopy<
              member_type, KokkosBatched::Trans::Transpose>::invoke(member, Pi,
                                                                    PTri);
          KokkosBatched::TeamGemm<
              member_type, KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, Mi, T,
                                                            0.0, MTi);
          member.team_barrier();
          KokkosBatched::TeamLU<
              member_type, KokkosBatched::Algo::LU::Unblocked>::invoke(member,
                                                                       MTi);
          member.team_barrier();
          KokkosBatched::TeamSolveLU<
              member_type, KokkosBatched::Trans::Transpose,
              KokkosBatched::Algo::SolveLU::Unblocked>::invoke(member, MTi,
                                                               PTri);
          member.team_barrier();
          // move the solution into the solution vector Di (needs to be
          // transposed)
          KokkosBatched::TeamCopy<
              member_type, KokkosBatched::Trans::Transpose>::invoke(member,
                                                                    PTri, Di);
        });
    MandelToVoigt<ExeSpace>(D);
    Kokkos::fence();
    return D;
  }

  template <typename ViewT, typename ProbeViewT>
  auto ComputeProbeU(ViewT U, ProbeViewT probe, ViewT probe_U, double h)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 3 &&
                              Kokkos::is_view<ProbeViewT>::value &&
                              ProbeViewT::rank == 2,
                          void>
  {
    static_assert(Kokkos::SpaceAccessibility<
                      typename ViewT::memory_space,
                      typename ProbeViewT::memory_space>::accessible,
                  "ViewT and ProbeViewT must have the same memory space");
    using exe_space = typename ViewT::execution_space;
    KOKKOS_ASSERT(U.extent(0) == probe_U.extent(0));
    KOKKOS_ASSERT(U.extent(1) == probe_U.extent(1));
    KOKKOS_ASSERT(U.extent(2) == probe_U.extent(2));
    KOKKOS_ASSERT(U.extent(1) == probe.extent(0));
    KOKKOS_ASSERT(U.extent(2) == probe.extent(1));

    Kokkos::parallel_for(
        "calc probing vec",
        Kokkos::MDRangePolicy<exe_space, Kokkos::Rank<3>,
                              Kokkos::IndexType<size_t>>(
            {0ul, 0ul, 0ul}, {U.extent(0), U.extent(1), U.extent(2)}),
        KOKKOS_LAMBDA(int j, int k, int l) {
          probe_U(j, k, l) = U(j, k, l) + h * probe(k, l);
        });
  }

  template <typename ViewT>
  auto ComputeF(ViewT R, ViewT U, ViewT F)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 3,
                          void>

  {
    using exe_space = typename ViewT::execution_space;
    KOKKOS_ASSERT((R.extent(0) == U.extent(0)) &&
                  (R.extent(1) == U.extent(1)) && (R.extent(2) == U.extent(2)));
    KOKKOS_ASSERT((R.extent(0) == F.extent(0)) &&
                  (R.extent(1) == F.extent(1)) && (R.extent(2) == F.extent(2)));

    auto team_policy =
        Kokkos::TeamPolicy<exe_space>(U.extent(0), Kokkos::AUTO());
    using member_type = typename decltype(team_policy)::member_type;
    Kokkos::parallel_for(
        "F=RU", team_policy, KOKKOS_LAMBDA(member_type member) {
          auto j = member.league_rank();
          auto R_j = Kokkos::subview(R, j, Kokkos::ALL(), Kokkos::ALL());
          auto U_j = Kokkos::subview(U, j, Kokkos::ALL(), Kokkos::ALL());
          auto F_j = Kokkos::subview(F, j, Kokkos::ALL(), Kokkos::ALL());
          KokkosBatched::TeamGemm<
              member_type, KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, R_j,
                                                            U_j, 0.0, F_j);
        });
  }

  // W: work array should be same size as F
  template <typename ViewT>
  auto ComputeFIncrement(ViewT F, ViewT Fnew, ViewT W, ViewT Fincrement)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 3,
                          void>
  {
    using exe_space = typename ViewT::execution_space;
    KOKKOS_ASSERT((F.extent(0) == Fnew.extent(0)) &&
                  (F.extent(1) == Fnew.extent(1)) &&
                  (F.extent(2) == Fnew.extent(2)));
    KOKKOS_ASSERT((F.extent(0) == Fincrement.extent(0)) &&
                  (F.extent(1) == Fincrement.extent(1)) &&
                  (F.extent(2) == Fincrement.extent(2)));
    KOKKOS_ASSERT((F.extent(0) == W.extent(0)) &&
                  (F.extent(1) == W.extent(1)) && (F.extent(2) == W.extent(2)));
    Kokkos::parallel_for(
        "Finc", Kokkos::RangePolicy<exe_space>(0, F.extent(0)),
        KOKKOS_LAMBDA(int j) {
          auto F_probe_j =
              Kokkos::subview(Fnew, j, Kokkos::ALL(), Kokkos::ALL());
          auto F_increment_j =
              Kokkos::subview(Fincrement, j, Kokkos::ALL(), Kokkos::ALL());
          auto F_j = Kokkos::subview(F, j, Kokkos::ALL(), Kokkos::ALL());
          // Reuse the Probing_U memory to store the inverse of F
          auto Finv_j = Kokkos::subview(W, j, Kokkos::ALL(), Kokkos::ALL());
          // compute the inverse of F
          Invert<3>(F_j, Finv_j);
          // compute the increment in F
          KokkosBatched::SerialGemm<
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Trans::NoTranspose,
              KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, F_probe_j,
                                                            Finv_j, 0.0,
                                                            F_increment_j);
        });
  }

  template <typename ExeSpace, typename ComputePK2Func>
  struct StressFiniteDifferenceFunc
  {
    using memory_space = typename ExeSpace::memory_space;

    explicit StressFiniteDifferenceFunc(
        ComputePK2Func compute_stress,
        Kokkos::View<Scalar * [6], memory_space> current_stress,
        double h = 1e-6)
        : compute_pk2_stress_(compute_stress), h_(h)
    {
      current_pk2_stress_ =
          Kokkos::View<Scalar * [6], typename ExeSpace::memory_space>(
              "current pk2 stress", current_stress.extent(0));
      Kokkos::deep_copy(current_pk2_stress_, current_stress);
    }

    auto operator()(Kokkos::View<Scalar * [3][3], memory_space> F,
                    Kokkos::View<Scalar * [3][3], memory_space> R,
                    Kokkos::View<Scalar * [3][3], memory_space> U) noexcept
        -> Kokkos::View<Scalar * [6][6], memory_space>
    {
      assert(F.extent(0) == R.extent(0));
      assert(R.extent(0) == U.extent(0));
      assert(R.extent(0) == current_pk2_stress_.extent(0));
      Kokkos::View<Scalar * [6][6], memory_space> dPK2dU("dPK2dU", R.extent(0));
      Kokkos::View<Scalar * [3][3], memory_space> probing_F("probe F",
                                                            R.extent(0));
      Kokkos::View<Scalar * [3][3], memory_space> probing_F_increment(
          "probe F increment", R.extent(0));
      Kokkos::View<Scalar * [3][3], memory_space> probing_U("probe F",
                                                            R.extent(0));
      auto probes = ComputeProbeVectors<ExeSpace>();
      assert(probes.extent(0) == 6);
      // probing the 6 directions
      for (int i = 0; i < 6; ++i)
      {
        // 1) Noting that we want to perturb U, but not the full F we compute
        // the the new value of U as U+probe
        auto probe = Kokkos::subview(probes, i, Kokkos::ALL(), Kokkos::ALL());
        ComputeProbeU(U, probe, probing_U, h_);
        // 2) compute Fprobe as R@(U+probe).
        ComputeF(R, probing_U, probing_F);
        // 3) compute the increment in F as Finc = Fprobe@F^-1
        // we do not use team parallelism here because current implementation of
        // Invert is only implemented for the first order parallelism
        // reuse probing_U as the work array for ComputeFIncrement
        ComputeFIncrement(F, probing_F, probing_U, probing_F_increment);
        auto updated_stress =
            compute_pk2_stress_(probing_F, probing_F_increment);
        Kokkos::fence();
        Kokkos::parallel_for(
            "finite difference",
            Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<2>,
                                  Kokkos::IndexType<size_t>>(
                {0ul, 0ul}, {U.extent(0), 6ul}),
            KOKKOS_LAMBDA(int j, int k) {
              dPK2dU(j, k, i) =
                  (updated_stress(j, k) - current_pk2_stress_(j, k)) / h_;
            });
      }
      VoigtToMandel<ExeSpace>(dPK2dU);
      return dPK2dU;
    }

    ComputePK2Func compute_pk2_stress_;
    Kokkos::View<Scalar * [6], typename ExeSpace::memory_space>
        current_pk2_stress_;
    double h_;
  };

  // deduction guide to test the execution space
  template <typename Func, typename View>
  StressFiniteDifferenceFunc(Func, View, double)
      -> StressFiniteDifferenceFunc<typename View::execution_space, Func>;

  template <typename ExeSpace, typename ComputePK2Func>
  struct StressCentralDifferenceFunc
  {
    using exe_space = ExeSpace;
    using memory_space = typename ExeSpace::memory_space;

    explicit StressCentralDifferenceFunc(ComputePK2Func compute_stress,
                                         double h = 1e-6)
        : compute_pk2_stress_(compute_stress), h_(h)
    {
    }

    auto operator()(Kokkos::View<Scalar * [3][3], memory_space> F,
                    Kokkos::View<Scalar * [3][3], memory_space> R,
                    Kokkos::View<Scalar * [3][3], memory_space> U) noexcept
        -> Kokkos::View<Scalar * [6][6], memory_space>
    {
      // allocate an array for the current stress if it hasn't already been
      // allocated
      if (current_pk2_stress_.extent(0) == 0)
      {
        current_pk2_stress_ =
            Kokkos::View<Scalar * [6], typename ExeSpace::memory_space>(
                "current pk2 stress", F.extent(0));
      }
      assert(F.extent(0) == R.extent(0));
      assert(R.extent(0) == U.extent(0));
      assert(R.extent(0) == current_pk2_stress_.extent(0));
      Kokkos::View<Scalar * [6][6], memory_space> dPK2dU("dPK2dU", R.extent(0));
      Kokkos::View<Scalar * [3][3], memory_space> probing_F("probe F",
                                                            R.extent(0));
      Kokkos::View<Scalar * [3][3], memory_space> probing_F_increment(
          "probe F increment", R.extent(0));
      Kokkos::View<Scalar * [3][3], memory_space> probing_U("probe F",
                                                            R.extent(0));
      auto probes = ComputeProbeVectors<ExeSpace>();
      assert(probes.extent(0) == 6);
      // probing the 6 directions
      for (int i = 0; i < 6; ++i)
      {
        auto probe = Kokkos::subview(probes, i, Kokkos::ALL(), Kokkos::ALL());
        ComputeProbeU(U, probe, probing_U, h_);
        ComputeF(R, probing_U, probing_F);
        ComputeFIncrement(F, probing_F, probing_U, probing_F_increment);
        Kokkos::deep_copy(current_pk2_stress_,
                          compute_pk2_stress_(probing_F, probing_F_increment));
        ComputeProbeU(U, probe, probing_U, -1 * h_);
        ComputeF(R, probing_U, probing_F);
        ComputeFIncrement(F, probing_F, probing_U, probing_F_increment);
        auto updated_stress_negative =
            compute_pk2_stress_(probing_F, probing_F_increment);
        Kokkos::fence();
        Kokkos::parallel_for(
            "finite difference",
            Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<2>,
                                  Kokkos::IndexType<size_t>>(
                {0ul, 0ul}, {U.extent(0), 6ul}),
            KOKKOS_LAMBDA(int j, int k) {
              dPK2dU(j, k, i) =
                  (current_pk2_stress_(j, k) - updated_stress_negative(j, k)) /
                  (2 * h_);
            });
      }
      VoigtToMandel<ExeSpace>(dPK2dU);
      return dPK2dU;
    }

    ComputePK2Func compute_pk2_stress_;
    Kokkos::View<Scalar * [6], typename ExeSpace::memory_space>
        current_pk2_stress_;
    double h_;
  };

  // deduction guide to test the execution space
  template <typename Func>
  StressCentralDifferenceFunc(Func, double)
      -> StressCentralDifferenceFunc<typename Func::exe_space, Func>;



  template <typename ExeSpace>
  void ConvertTLStiffnessToULStiffness(
      Kokkos::View<Scalar * [3][3], typename ExeSpace::memory_space> F,
      Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space> D)
  {
    assert(F.extent(0) == D.extent(0));
    Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space> C(
        "UL material stiffness", F.extent(0));
    Kokkos::deep_copy(C, 0.0);
    Kokkos::parallel_for(
        "TL2UL", Kokkos::RangePolicy<ExeSpace>(0, F.extent(0)),
        KOKKOS_LAMBDA(int i) {
          auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
          auto Ci = Kokkos::subview(C, i, Kokkos::ALL(), Kokkos::ALL());
          auto Di = Kokkos::subview(D, i, Kokkos::ALL(), Kokkos::ALL());
          for (int m = 0; m < 3; ++m)
          {
            for (int n = m; n < 3; ++n)
            {
              for (int p = 0; p < 3; ++p)
              {
                for (int q = p; q < 3; ++q)
                {
                  auto idx_mn = TensorIndex2VoigtIndex(m, n);
                  auto idx_pq = TensorIndex2VoigtIndex(p, q);
                  double total = 0.0;
                  for (int j = 0; j < 3; ++j)
                  {
                    for (int k = 0; k < 3; ++k)
                    {
                      for (int r = 0; r < 3; ++r)
                      {
                        for (int s = 0; s < 3; ++s)
                        {
                          total += Fi(m, j) * Fi(n, k) * Fi(p, r) * Fi(q, s) *
                              Di(TensorIndex2VoigtIndex(j, k),
                                 TensorIndex2VoigtIndex(r, s));
                        }
                      }
                    }
                  }
                  Ci(idx_mn, idx_pq) += total;
                }
              }
            }
          }
          KokkosBatched::SerialScale::invoke(1.0 / Determinant<3>(Fi), Ci);
        });
    Kokkos::deep_copy(D, C);
  }

}  // namespace mumfim

#endif  // MUMFIM_SRC_MUMFIM_MICROSCALE_MATERIALSTIFFNESS_H
