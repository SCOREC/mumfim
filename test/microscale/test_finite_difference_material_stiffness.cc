#include <mumfim/microscale/MicroTypeDefinitions.h>

#include <Catch2/catch.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
// #include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemm_Team_Impl.hpp>
// #include <KokkosBatched_InverseLU_Decl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
// #include <KokkosBatched_LU_Serial_Impl.hpp>
#include <mumfim/microscale/BatchedNeohookeanAnalysis.h>
#include <mumfim/microscale/BatchedRVEAnalysis.h>
#include <mumfim/microscale/PolarDecomposition.h>

#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Copy_Impl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
#include <KokkosBatched_Scale_Decl.hpp>
#include <KokkosBatched_Scale_Impl.hpp>
#include <KokkosBatched_SolveLU_Decl.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <chrono>

namespace mumfim
{
  /**
   * @brief Compute the determinant of a 2x2 or 3x3 matrix
   * @param[in] F the matrix
   * @return the determinant
   * @note this function is KOKKOS_INLINE and can be used on the device
   */
  template <int dim, typename ViewT>
  KOKKOS_INLINE_FUNCTION constexpr auto Determinant(ViewT F)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 2,
                          double>
  {
    static_assert(dim == 2 || dim == 3, "only 2d and 3d implemented");
    assert(F.extent(0) == dim && F.extent(1) == dim);
    if constexpr (dim == 2)
    {
      return (F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0));
    }
    else if constexpr (dim == 3)
    {
      // return (F(0, 0) * (F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1)) -
      //         F(0, 1) * (F(1, 0) * F(2, 2) - F(1, 2) * F(2, 0)) +
      //         F(0, 2) * (F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0)));
      return (F(0, 0) * F(1, 1) * F(2, 2) + F(0, 1) * F(1, 2) * F(2, 0) +
              F(0, 2) * F(1, 0) * F(2, 1) -
              (F(0, 2) * F(1, 1) * F(2, 0) + F(0, 1) * F(1, 0) * F(2, 2) +
               F(0, 0) * F(1, 2) * F(2, 1)));
    }
  }

  /**
   * @brief Compute the inverse of a 2x2 or 3x3 matrix
   * @param[in] F the matrix
   * @param[out] F_inv the inverse of F
   * @return the determinant of F
   * @note this function is KOKKOS_INLINE and can be used on the device
   */
  template <int dim, typename ViewT>
  KOKKOS_INLINE_FUNCTION constexpr auto Invert(ViewT F, ViewT invF)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 2,
                          void>
  {
    static_assert(dim == 2 || dim == 3, "only 2d and 3d implemented");
    assert(F.extent(0) == dim && F.extent(1) == dim);
    if constexpr (dim == 2)
    {
      double detF = Determinant<2>(F);
      KOKKOS_ASSERT(detF > 0.0);
      invF(0, 0) = F(1, 1) / detF;
      invF(0, 1) = -F(0, 1) / detF;
      invF(1, 0) = -F(1, 0) / detF;
      invF(1, 1) = F(0, 0) / detF;
    }
    else if constexpr (dim == 3)
    {
      double detF = Determinant<3>(F);
      KOKKOS_ASSERT(detF > 0.0);
      invF(0, 0) = (F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1)) / detF;
      invF(0, 1) = (F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2)) / detF;
      invF(0, 2) = (F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1)) / detF;
      invF(1, 0) = (F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2)) / detF;
      invF(1, 1) = (F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0)) / detF;
      invF(1, 2) = (F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2)) / detF;
      invF(2, 0) = (F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0)) / detF;
      invF(2, 1) = (F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1)) / detF;
      invF(2, 2) = (F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0)) / detF;
    }
  }

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
  auto ComputeProbeVectors(double h)
      -> Kokkos::View<Scalar[6][3][3], MemorySpace>
  {
    Kokkos::View<Scalar[6][3][3], MemorySpace> T("T");
    auto T_h = Kokkos::create_mirror_view(Kokkos::HostSpace{}, T);
    Kokkos::deep_copy(T_h, 0.0);
    T_h(0, 0, 0) = h;
    T_h(1, 1, 1) = h;
    T_h(2, 2, 2) = h;
    T_h(3, 1, 1) = h;
    T_h(3, 1, 2) = h;
    T_h(3, 2, 1) = h;
    T_h(3, 2, 2) = h;
    T_h(4, 0, 0) = h;
    T_h(4, 0, 2) = h;
    T_h(4, 2, 0) = h;
    T_h(4, 2, 2) = h;
    T_h(5, 0, 0) = h;
    T_h(5, 0, 1) = h;
    T_h(5, 1, 0) = h;
    T_h(5, 1, 1) = h;
    Kokkos::deep_copy(T, T_h);
    return T;
  }

  template <typename ExeSpace>
  void MandelToVoigt(
      Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space> matrix)
  {
    auto sr2 = sqrt(2);
    Kokkos::parallel_for(
        "mandel to voigt",
        Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>>(
            {0, 0, 0}, {matrix.extent(0), matrix.extent(1), matrix.extent(2)}),
        KOKKOS_LAMBDA(int i, int j, int k) {
          // lower right
          if (j > 2 && k > 2)
          {
            matrix(i, j, k) = matrix(i, j, k) / 2;
          }
          // upper right
          else if (j < 3 && k > 2)
          {
            matrix(i, j, k) = matrix(i, j, k) / sr2;
          }
          // lower left
          else if (j > 2)
          {
            matrix(i, j, k) = matrix(i, j, k) / sr2;
          }
        });
  }

  template <typename ExeSpace>
  void VoigtToMandel(
      Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space> matrix)
  {
    auto sr2 = sqrt(2);
    Kokkos::parallel_for(
        "voigt to mandel",
        Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>,
                              Kokkos::IndexType<size_t>>(
            {0ul, 0ul, 0ul},
            {matrix.extent(0), matrix.extent(1), matrix.extent(2)}),
        KOKKOS_LAMBDA(int i, int j, int k) {
          // lower right
          if (j > 2 && k > 2)
          {
            matrix(i, j, k) = 2 * matrix(i, j, k);
          }
          // upper right
          else if (j < 3 && k > 2)
          {
            matrix(i, j, k) = sr2 * matrix(i, j, k);
          }
          // lower left
          else if (j > 2)
          {
            matrix(i, j, k) = sr2 * matrix(i, j, k);
          }
        });
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

  template <typename ExeSpace>
  auto TLtoULStiffness(
      Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space> dPK2dE)
      -> Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space>
  {
    Kokkos::parallel_for(
        "TL2UL stiffness", Kokkos::RangePolicy<ExeSpace>(0, dPK2dE.extent(0)),
        KOKKOS_LAMBDA(size_t idx, size_t j, size_t k, size_t l, size_t m){

        });
  }

  template <typename ExeSpace, typename ComputePK2Func>
  struct StressFiniteDifferenceFunc
  {
    using exe_space = ExeSpace;
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
      auto probes = ComputeProbeVectors<ExeSpace>(h_);
      assert(probes.extent(0) == 6);
      // probing the 6 directions
      for (int i = 0; i < 6; ++i)
      {
        // 1) Noting that we want to perturb U, but not the full F we compute
        // the the new value of U as U+probe
        Kokkos::parallel_for(
            "calc probing vec",
            Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>,
                                  Kokkos::IndexType<size_t>>(
                {0ul, 0ul, 0ul}, {R.extent(0), 3ul, 3ul}),
            KOKKOS_LAMBDA(int j, int k, int l) {
              probing_U(j, k, l) = U(j, k, l) + probes(i, k, l);
            });
        // 2) compute Fprobe as R@(U+probe).
        auto team_policy =
            Kokkos::TeamPolicy<ExeSpace>(U.extent(0), Kokkos::AUTO());
        using member_type = typename decltype(team_policy)::member_type;
        Kokkos::parallel_for(
            "F=RU", team_policy, KOKKOS_LAMBDA(member_type member) {
              auto j = member.league_rank();
              auto R_j = Kokkos::subview(R, j, Kokkos::ALL(), Kokkos::ALL());
              auto U_j =
                  Kokkos::subview(probing_U, j, Kokkos::ALL(), Kokkos::ALL());
              auto F_j =
                  Kokkos::subview(probing_F, j, Kokkos::ALL(), Kokkos::ALL());
              KokkosBatched::TeamGemm<
                  member_type, KokkosBatched::Trans::NoTranspose,
                  KokkosBatched::Trans::NoTranspose,
                  KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0,
                                                                R_j, U_j, 0.0,
                                                                F_j);
            });
        // 3) compute the increment in F as Finc = Fprobe@F^-1
        // we do not use team parallelism here because current implementation of
        // Invert is only implemented for the first order parallelism
        Kokkos::parallel_for(
            "Finc", Kokkos::RangePolicy<ExeSpace>(0, U.extent(0)),
            KOKKOS_LAMBDA(int j) {
              auto F_probe_j =
                  Kokkos::subview(probing_F, j, Kokkos::ALL(), Kokkos::ALL());
              auto F_increment_j = Kokkos::subview(
                  probing_F_increment, j, Kokkos::ALL(), Kokkos::ALL());
              auto F_j = Kokkos::subview(F, j, Kokkos::ALL(), Kokkos::ALL());
              // Reuse the Probing_U memory to store the inverse of F
              auto Finv_j =
                  Kokkos::subview(probing_U, j, Kokkos::ALL(), Kokkos::ALL());
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

  /**
   * @brief Convert a stress tensor from Voigt to matrix form
   * @param[in] stress_voigt the stress tensor in Voigt form
   * @param[out] stress_matrix the stress tensor in matrix form
   * @note this function is KOKKOS_INLINE and can be used on the device
   */
  template <typename ViewT, typename ViewT2>
  KOKKOS_INLINE_FUNCTION constexpr auto VoigtToMatrix(ViewT stress_voigt,
                                                      ViewT2 stress_matrix)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 1 &&
                              Kokkos::is_view<ViewT2>::value &&
                              ViewT2::rank == 2,
                          void>
  {
    KOKKOS_ASSERT(stress_voigt.extent(0) == 6);
    KOKKOS_ASSERT(stress_matrix.extent(0) == 3 && stress_matrix.extent(1) == 3);
    stress_matrix(0, 0) = stress_voigt(0);
    stress_matrix(1, 1) = stress_voigt(1);
    stress_matrix(2, 2) = stress_voigt(2);
    stress_matrix(2, 1) = stress_voigt(3);
    stress_matrix(1, 2) = stress_voigt(3);
    stress_matrix(0, 2) = stress_voigt(4);
    stress_matrix(2, 0) = stress_voigt(4);
    stress_matrix(0, 1) = stress_voigt(5);
    stress_matrix(1, 0) = stress_voigt(5);
  }

  /**
   * @brief Convert a stress tensor from matrix to Voigt form
   * @param[in] stress_matrix the stress tensor in matrix form
   * @param[out] stress_voigt the stress tensor in Voigt form
   * @note this function is KOKKOS_INLINE and can be used on the device
   */
  template <typename ViewT, typename ViewT2>
  KOKKOS_INLINE_FUNCTION constexpr auto MatrixToVoigt(ViewT stress_matrix,
                                                      ViewT2 stress_voigt)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 2 &&
                              Kokkos::is_view<ViewT2>::value &&
                              ViewT2::rank == 1,
                          void>
  {
    stress_voigt(0) = stress_matrix(0, 0);
    stress_voigt(1) = stress_matrix(1, 1);
    stress_voigt(2) = stress_matrix(2, 2);
    stress_voigt(3) = stress_matrix(2, 1);
    stress_voigt(4) = stress_matrix(0, 2);
    stress_voigt(5) = stress_matrix(1, 0);
  }

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

  KOKKOS_INLINE_FUNCTION
  size_t TensorIndex2VoigtIndex(size_t i, size_t j)
  {
    // 11, 22, 33, 23, 13, 12
    if (i == j)
    {
      return i;
    }
    else if ((i == 1 && j == 2) || (i == 2 && j == 1))
    {
      return 3;
    }
    else if ((i == 0 && j == 2) || (i == 2 && j == 0))
    {
      return 4;
    }
    else if ((i == 0 && j == 1) || (i == 1 && j == 0))
    {
      return 5;
    }
    else
    {
      KOKKOS_ASSERT(i < 3 && j < 3);
    }
    return -1;
  }

  KOKKOS_INLINE_FUNCTION
  size_t TensorIndex2VoigtIndexNoTranspose(size_t i, size_t j)
  {
    // 11, 22, 33, 23, 13, 12
    if (i == j)
    {
      return i;
    }
    else if ((i == 1 && j == 2))
    {
      return 3;
    }
    else if ((i == 0 && j == 2))
    {
      return 4;
    }
    else if ((i == 0 && j == 1))
    {
      return 5;
    }
    return 99999;
  }

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

  template <typename Analysis>
  struct BatchedAnalysisGetPK2StressFunc
  {
    using memory_space = typename Analysis::memory_space;
    using exe_space = typename Analysis::exe_space;

    explicit BatchedAnalysisGetPK2StressFunc(Analysis & analysis)
        : analysis_(analysis)
        , F_increment_("deformation grad increment", analysis_.GetNumRVEs())
        , stress_("stress", analysis_.GetNumRVEs())
    {
    }

    // That way we can work with codes that use either incrementalor total
    // approach
    auto operator()(Kokkos::View<Scalar * [3][3], memory_space> F,
                    Kokkos::View<Scalar * [3][3], memory_space> F_increment,
                    bool update_coords = false) noexcept
        -> Kokkos::View<Scalar * [6], memory_space>
    {
      assert(F_increment.extent(0) == analysis_.GetNumRVEs());
      assert(F.extent(0) == analysis_.GetNumRVEs());
      // FIXME pass in here Fincrement
      Kokkos::deep_copy(F_increment_.d_view, F_increment);
      F_increment_.template modify<exe_space>();
      analysis_.run(F_increment_, stress_, update_coords);
      stress_.template sync<exe_space>();
      // Batched analysis computes the cauchy stress, we need to convert this to
      // the PK2 stress
      ConvertCauchyToPK2<exe_space>(F, stress_.d_view);
      return stress_.d_view;
    }

    Analysis & analysis_;
    Kokkos::DualView<Scalar * [3][3], memory_space> F_increment_;
    Kokkos::DualView<Scalar * [6], memory_space> stress_;
  };
}  // namespace mumfim

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
  constexpr int num_tests = 1;
  Kokkos::View<double * [3][3], Kokkos::HostSpace> F("F", num_tests);
  Kokkos::View<double * [6], Kokkos::HostSpace> PK2("PK2", num_tests);
  Kokkos::fill_random(F, random, 1.0);
  Kokkos::parallel_for(
      F.extent(0), KOKKOS_LAMBDA(int i) {
        F(i, 0, 0) += 1.0;
        F(i, 1, 1) += 1.0;
        F(i, 2, 2) += 1.0;
      });
  Kokkos::parallel_for(
      F.extent(0), KOKKOS_LAMBDA(int i) {
        auto Fi = Kokkos::subview(F, i, Kokkos::ALL(), Kokkos::ALL());
        KOKKOS_ASSERT(mumfim::Determinant<3>(Fi) > 0.0);
      });

  auto t1 = std::chrono::steady_clock::now();
  // mumfim::computeDPK2dE<Kokkos::Serial>(F, dPK2dU);
  auto t2 = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration<double>(t2 - t1).count();
  std::cout << "Time taken by function: " << duration << " seconds"
            << std::endl;
  mumfim::BatchedNeohookeanAnalysis analysis(F.extent(0), 1E3, 0.3);
  using exe_space = typename decltype(analysis)::exe_space;
  mumfim::BatchedAnalysisGetPK2StressFunc compute_pk2_stress{analysis};
  //  LinearStress<exe_space> compute_pk2_stress{};
  // CubicStress<exe_space> compute_pk2_stress{};
  auto [U, R] = mumfim::PolarDecomposition<exe_space>(F);
  // apply the full deformation gradient as the increment
  // in the first timestep
  auto stress = compute_pk2_stress(F, F, true);
  for (size_t i = 0; i < 3; ++i)
  {
    std::cout << "(";
    for (size_t j = 0; j < 3; ++j)
    {
      std::cout << F(0, i, j);
      if (j < 2) std::cout << ",";
    }
    std::cout << "),";
  }

  std::cout << "\n";
  mumfim::StressFiniteDifferenceFunc finite_difference{compute_pk2_stress,
                                                       stress, 1E-7};
  // auto result = finite_difference(F, R, U);
  auto result = mumfim::ComputeDPK2dE<exe_space>(F, finite_difference);
  for (size_t i = 0; i < result.extent(1); ++i)
  {
    for (size_t j = 0; j < result.extent(2); ++j)
    {
      auto r = result(0, i, j) < 1E-5 ? 0.0 : result(0, i, j);
      std::cout << r << " ";
    }
    std::cout << "\n";
  }
  std::cout << "---------------------\n";
  mumfim::ConvertTLStiffnessToULStiffness<exe_space>(F, result);
  for (size_t i = 0; i < result.extent(1); ++i)
  {
    for (size_t j = 0; j < result.extent(2); ++j)
    {
      auto r = result(0, i, j) < 1E-5 ? 0.0 : result(0, i, j);
      std::cout << r << " ";
    }
    std::cout << "\n";
  }
  Kokkos::DualView<double * [6][6], typename exe_space::memory_space> stiffness(
      "stiffness", F.extent(0));
  analysis.computeMaterialStiffness(stiffness);
  stiffness.sync_host();
  auto stiffness_h = stiffness.h_view;
  std::cout << "---------------------\n";
  for (size_t i = 0; i < stiffness_h.extent(1); ++i)
  {
    for (size_t j = 0; j < stiffness_h.extent(2); ++j)
    {
      auto r = stiffness_h(0, i, j) < 1E-5 ? 0.0 : stiffness_h(0, i, j);
      std::cout << r << " ";
    }
    std::cout << "\n";
  }
}
