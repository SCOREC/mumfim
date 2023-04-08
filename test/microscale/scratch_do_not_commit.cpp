template <typename ExeSpace, typename Func>
auto StressFiniteDifference(
    Kokkos::View<Scalar * [3][3], typename ExeSpace::memory_space> R,
    Kokkos::View<Scalar * [3][3], typename ExeSpace::memory_space> U,
    Kokkos::View<Scalar * [6], typename ExeSpace::memory_space>
    current_pk2_stress,
    Func compute_pk2_stress,
    double h = 1E-8)
-> Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space>
{
  using memory_space = typename ExeSpace::memory_space;
  assert(R.extent(0) == U.extent(0));
  assert(R.extent(0) == current_pk2_stress.extent(0));
  Kokkos::View<Scalar * [6][6], memory_space> dPK2dU("dPK2dU", R.extent(0));
  Kokkos::View<Scalar * [3][3], memory_space> probing_F("probe F",
                                                        R.extent(0));
  Kokkos::View<Scalar * [3][3], memory_space> probing_U("probe F",
                                                        R.extent(0));
  auto probes = ComputeProbeVectors<ExeSpace>();
  assert(probes.extent(0) == 6);
  // probing the 6 directions
  for (int i = 0; i < 6; ++i)
  {
    Kokkos::parallel_for(
        "calc probing vec",
        Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>,
        Kokkos::IndexType<size_t>>(
            {0ul, 0ul, 0ul}, {R.extent(0), 3ul, 3ul}),
        KOKKOS_LAMBDA(int j, int k, int l) {
      probing_U(j, k, l) = U(j, k, l) + h * probes(i, k, l);
    });
    auto team_policy =
        Kokkos::TeamPolicy<ExeSpace>(U.extent(0), Kokkos::AUTO());
    using member_type = typename decltype(team_policy)::member_type;
    Kokkos::parallel_for(
        "F=RU", team_policy, KOKKOS_LAMBDA(member_type member) {
      auto i = member.league_rank();
      auto Ri = Kokkos::subview(R, i, Kokkos::ALL(), Kokkos::ALL());
      auto Ui =
          Kokkos::subview(probing_U, i, Kokkos::ALL(), Kokkos::ALL());
      auto Fi =
          Kokkos::subview(probing_F, i, Kokkos::ALL(), Kokkos::ALL());
      KokkosBatched::TeamGemm<
          member_type, KokkosBatched::Trans::NoTranspose,
          KokkosBatched::Trans::NoTranspose,
          KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, Ri,
                                                        Ui, 0.0, Fi);
    });
    auto updated_pk2_stress = compute_pk2_stress(probing_F);
    Kokkos::fence();
    Kokkos::parallel_for(
        "finite difference",
        Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<2>,
        Kokkos::IndexType<size_t>>({0ul, 0ul},
        {U.extent(0), 6ul}),
        KOKKOS_LAMBDA(int j, int k) {
      dPK2dU(j, k, i) =
          (updated_pk2_stress(j, k) - current_pk2_stress(j, k)) / h;
    });
  }
  return dPK2dU;
}
