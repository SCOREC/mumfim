#ifndef MUMFIM_BATCHED_COMPUTE_ORIENTATION_TENSOR_H
#define MUMFIM_BATCHED_COMPUTE_ORIENTATION_TENSOR_H
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include "bioPackedData.h"
namespace mumfim
{
  template <typename Scalar, typename Ordinal, typename CoordinatesType, typename ConnectivityType, typename ExeSpace>
  struct OrientationTensor
  {
    struct TagCompute3D
    {
    };
    struct TagCompute2D
    {
    };
    //using CoordinatesType = PackedData<Scalar*, ExeSpace>;
    //using ConnectivityType = PackedData<Ordinal*[2], ExeSpace>;

    using OrientationType = Kokkos::DualView<Scalar * [3][3], ExeSpace>;
    using NormalType = Kokkos::DualView<Scalar * [3], ExeSpace>;
    using OuterPolicyType3D = Kokkos::TeamPolicy<ExeSpace, TagCompute3D>;
    using OuterPolicyType2D = Kokkos::TeamPolicy<ExeSpace, TagCompute2D>;
    using member_type_3D = typename OuterPolicyType3D::member_type;
    using member_type_2D = typename OuterPolicyType2D::member_type;
    using OrientationDeviceViewType =
        Kokkos::View<typename OrientationType::traits::data_type,
                     typename OrientationType::traits::array_layout,
                     ExeSpace>;
    using NormalDeviceViewType =
        Kokkos::View<typename NormalType::traits::data_type,
                     typename NormalType::traits::array_layout,
                     ExeSpace>;
    using OrientationScratchPad =
        Kokkos::View<Scalar[3][3],
                     typename ExeSpace::scratch_memory_space,
                     Kokkos::MemoryUnmanaged>;
    using NormalScratchPad =
        Kokkos::View<Scalar[3],
                     typename ExeSpace::scratch_memory_space,
                     Kokkos::MemoryUnmanaged>;
    ConnectivityType connectivity_;
    CoordinatesType coordinates_;
    // Ordinal team_size_;
    typename OuterPolicyType3D::index_type team_size_;
    OrientationDeviceViewType omega_d_;
    NormalDeviceViewType normal_d_;
    OrientationTensor()= default;;
    OrientationTensor(ConnectivityType connectivity, CoordinatesType coordinates, Ordinal team_size)
        : connectivity_(connectivity)
        , coordinates_(coordinates)
        , team_size_(team_size)
    {
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagCompute2D, const member_type_2D & team_member) const
    {
      auto coordinates_row =
          coordinates_.template getRow<ExeSpace>(team_member.league_rank());
      auto connectivity_row =
          connectivity_.template getRow<ExeSpace>(team_member.league_rank());
      auto num_edges = connectivity_row.extent(0);
      OrientationScratchPad scratch_orientation(team_member.team_scratch(0));
      NormalScratchPad scratch_normal(team_member.team_scratch(0));
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, 0, 3),
                           [=](const int i) {
                             scratch_orientation(i, 0) = 0;
                             scratch_orientation(i, 1) = 0;
                             scratch_orientation(i, 2) = 0;
                           });
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team_member, 0, 3), [=](const int i) {
            scratch_normal(i) = normal_d_(team_member.league_rank(), i);
          });
      team_member.team_barrier();
      // normalize the scratch normal currently we do in serial
      Kokkos::single(Kokkos::PerTeam(team_member), [=]() {
        auto n = scratch_normal(0) * scratch_normal(0) +
                 scratch_normal(1) * scratch_normal(1) +
                 scratch_normal(2) * scratch_normal(2);
        scratch_normal(0) = scratch_normal(0) / n;
        scratch_normal(1) = scratch_normal(1) / n;
        scratch_normal(2) = scratch_normal(2) / n;
      });
      team_member.team_barrier();
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team_member, 0, num_edges), [=](const int i) {
            auto n1 = connectivity_row(i,0);
            auto n2 = connectivity_row(i,1);
            auto lx = coordinates_row(n2,0) - coordinates_row(n1,0);
            auto ly = coordinates_row(n2,1) - coordinates_row(n1,1);
            auto lz = coordinates_row(n2,2) - coordinates_row(n1,2);
            auto l = sqrt(lx * lx + ly * ly + lz * lz);
            // normalize the lx,ly,lz array
            lx = lx / l;
            ly = ly / l;
            lz = lz / l;
            // project the edge unit into plane defined by the normal
            auto dotn = lx * scratch_normal(0) + ly * scratch_normal(1) +
                        lz * scratch_normal(2);
            lx = lx - scratch_normal(0) * dotn;
            ly = ly - scratch_normal(1) * dotn;
            lz = lz - scratch_normal(2) * dotn;
            l = sqrt(lx * lx + ly * ly + lz * lz);
            lx = lx / l;
            ly = ly / l;
            lz = lz / l;
            Kokkos::atomic_add(&scratch_orientation(0, 0), lx * lx);
            Kokkos::atomic_add(&scratch_orientation(0, 1), lx * ly);
            Kokkos::atomic_add(&scratch_orientation(0, 2), lx * lz);
            Kokkos::atomic_add(&scratch_orientation(1, 0), ly * lx);
            Kokkos::atomic_add(&scratch_orientation(1, 1), ly * ly);
            Kokkos::atomic_add(&scratch_orientation(1, 2), ly * lz);
            Kokkos::atomic_add(&scratch_orientation(2, 0), lz * lx);
            Kokkos::atomic_add(&scratch_orientation(2, 1), lz * ly);
            Kokkos::atomic_add(&scratch_orientation(2, 2), lz * lz);
          });
      team_member.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, 0, 3),
                           [=](const int i) {
                             omega_d_(team_member.league_rank(), i, 0) =
                                 scratch_orientation(i, 0) / num_edges;
                             omega_d_(team_member.league_rank(), i, 1) =
                                 scratch_orientation(i, 1) / num_edges;
                             omega_d_(team_member.league_rank(), i, 2) =
                                 scratch_orientation(i, 2) / num_edges;
                           });
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagCompute3D, const member_type_3D & team_member) const
    {
      auto coordinates_row =
          coordinates_.template getRow<ExeSpace>(team_member.league_rank());
      auto connectivity_row =
          connectivity_.template getRow<ExeSpace>(team_member.league_rank());
      Ordinal num_edges = connectivity_row.extent(0);
      OrientationScratchPad scratch_orientation(team_member.team_scratch(0));
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, 0, 3),
                           [=](const int i) {
                             scratch_orientation(i, 0) = 0;
                             scratch_orientation(i, 1) = 0;
                             scratch_orientation(i, 2) = 0;
                           });
      team_member.team_barrier();
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team_member, 0, num_edges), [=](const int i) {
            auto n1 = connectivity_row(i,0);
            auto n2 = connectivity_row(i,1);
            auto lx = coordinates_row(n2,0) - coordinates_row(n1,0);
            auto ly = coordinates_row(n2,1) - coordinates_row(n1,1);
            auto lz = coordinates_row(n2,2) - coordinates_row(n1,2);
            auto l = sqrt(lx * lx + ly * ly + lz * lz);
            lx = lx / l;
            ly = ly / l;
            lz = lz / l;
            Kokkos::atomic_add(&scratch_orientation(0, 0), lx * lx);
            Kokkos::atomic_add(&scratch_orientation(0, 1), lx * ly);
            Kokkos::atomic_add(&scratch_orientation(0, 2), lx * lz);
            Kokkos::atomic_add(&scratch_orientation(1, 0), ly * lx);
            Kokkos::atomic_add(&scratch_orientation(1, 1), ly * ly);
            Kokkos::atomic_add(&scratch_orientation(1, 2), ly * lz);
            Kokkos::atomic_add(&scratch_orientation(2, 0), lz * lx);
            Kokkos::atomic_add(&scratch_orientation(2, 1), lz * ly);
            Kokkos::atomic_add(&scratch_orientation(2, 2), lz * lz);
          });
      team_member.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, 0, 3),
                           [=](const int i) {
                             omega_d_(team_member.league_rank(), i, 0) =
                                 scratch_orientation(i, 0) / num_edges;
                             omega_d_(team_member.league_rank(), i, 1) =
                                 scratch_orientation(i, 1) / num_edges;
                             omega_d_(team_member.league_rank(), i, 2) =
                                 scratch_orientation(i, 2) / num_edges;
                           });
    }
    void compute3D(OrientationType omega)
    {
      connectivity_.template sync<ExeSpace>();
      coordinates_.template sync<ExeSpace>();
      omega.template modify<ExeSpace>();
      // Kokkos::deep_copy(omega_d_, 0);
      omega_d_ = omega.template view<ExeSpace>();
      auto shared_bytes = OrientationScratchPad::shmem_size();
      Kokkos::parallel_for(
          OuterPolicyType3D(omega_d_.extent(0), team_size_)
              .set_scratch_size(0, Kokkos::PerTeam(shared_bytes)),
          *this);
      omega_d_ = OrientationDeviceViewType{};
    }
    void compute2D(NormalType normal, OrientationType omega)
    {
      connectivity_.template sync<ExeSpace>();
      coordinates_.template sync<ExeSpace>();
      normal.template sync<ExeSpace>();
      omega.template modify<ExeSpace>();
      omega_d_ = omega.template view<ExeSpace>();
      normal_d_ = normal.template view<ExeSpace>();
      Kokkos::deep_copy(omega_d_, 0);
      auto shared_bytes =
          OrientationScratchPad::shmem_size() + NormalScratchPad::shmem_size();
      Kokkos::parallel_for(
          OuterPolicyType2D(omega_d_.extent(0), team_size_)
              .set_scratch_size(0, Kokkos::PerTeam(shared_bytes)),
          *this);
      omega_d_ = OrientationDeviceViewType{};
      normal_d_ = NormalDeviceViewType{};
    }
  };
}  // namespace mumfim
#endif
