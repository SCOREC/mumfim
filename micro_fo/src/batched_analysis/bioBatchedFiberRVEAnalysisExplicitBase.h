#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_BASE_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_BASE_H__
#include <Kokkos_Core.hpp>
#include <array>
#include <limits>
#include "bioPackedData.h"
#include "bioRVE.h"
//#define TEAM_SIZE 512
#define TEAM_SIZE Kokkos::AUTO
namespace bio
{
  template <typename ExeSpace,
            typename Ordinal,
            typename T,
            typename T2,
            typename T3,
            int DIM = 3>
  void applyIncrementalDeformationToDisplacement(Ordinal number_rves,
                                                 T deformation_gradient,
                                                 T3 coordinates,
                                                 T2 displacements)
  {
    deformation_gradient.template sync<ExeSpace>();
    auto deformation_gradient_d =
        deformation_gradient.template view<ExeSpace>();
    coordinates.template sync<ExeSpace>();
    displacements.template sync<ExeSpace>();
    displacements.template modify<ExeSpace>();
    using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace>;
    using member_type = typename OuterPolicyType::member_type;
    OuterPolicyType outer_policy(number_rves, TEAM_SIZE);
    Kokkos::parallel_for(
        "applyIncrementalDfm", outer_policy,
        KOKKOS_LAMBDA(member_type team_member) {
          int rve = team_member.league_rank();
          auto displacements_row = displacements.template getRow<ExeSpace>(rve);
          auto coordinates_row = coordinates.template getRow<ExeSpace>(rve);
          auto deformation_idx = deformation_gradient_d.extent(0) > 1 ? rve : 0;
          auto deformation_gradient_row =
              Kokkos::subview(deformation_gradient_d, deformation_idx,
                              Kokkos::ALL(), Kokkos::ALL());
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team_member,
                                      displacements_row.extent(0) / DIM),
              [&](int i) {
                for (int k = 0; k < 3; ++k)
                {
                  displacements_row(i * DIM + k) =
                      (deformation_gradient_row(k, 0) - (k == 0 ? 1 : 0)) *
                          coordinates_row(i * DIM) +
                      (deformation_gradient_row(k, 1) - (k == 1 ? 1 : 0)) *
                          coordinates_row(i * DIM + 1) +
                      (deformation_gradient_row(k, 2) - (k == 2 ? 1 : 0)) *
                          coordinates_row(i * DIM + 2) +
                      displacements_row(i * DIM + k);
                }
              });
        });
  }
  template <typename ExeSpace,
            typename T,
            typename T2,
            typename T3,
            typename T4>
  void computeCauchyStress(T2 boundary_dofs,
                           T3 coordinates,
                           T3 force,
                           T stress,
                           T4 volume,
                           T4 scale_factor)
  {
    using Scalar = typename T::value_type;
    boundary_dofs.template sync<ExeSpace>();
    coordinates.template sync<ExeSpace>();
    force.template sync<ExeSpace>();
    volume.template sync<ExeSpace>();
    scale_factor.template sync<ExeSpace>();
    auto stress_d = stress.template view<ExeSpace>();
    auto volume_d = volume.template view<ExeSpace>();
    auto scale_factor_d = scale_factor.template view<ExeSpace>();
    // zero the trss vector because we will sum the values into the streses
    Kokkos::deep_copy(stress_d, 0);
    using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace>;
    using member_type = typename OuterPolicyType::member_type;
    OuterPolicyType outer_policy(stress.extent(0), TEAM_SIZE);
    Kokkos::parallel_for(
        "compute cauchy stress", outer_policy,
        KOKKOS_LAMBDA(member_type team_member) {
          int rve = team_member.league_rank();
          auto stress_row = Kokkos::subview(stress_d, rve, Kokkos::ALL());
          auto boundary_dof_row = boundary_dofs.template getRow<ExeSpace>(rve);
          auto coordinates_row = coordinates.template getRow<ExeSpace>(rve);
          auto force_row = force.template getRow<ExeSpace>(rve);
          auto loop_size = boundary_dof_row.extent(0) / 3;
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team_member, loop_size),
              [=](const int i) {
                auto dof1 = boundary_dof_row(i * 3);
                auto dof2 = boundary_dof_row(i * 3 + 1);
                auto dof3 = boundary_dof_row(i * 3 + 2);
                auto f1 = force_row(dof1);
                auto f2 = force_row(dof2);
                auto f3 = force_row(dof3);
                auto crd1 = coordinates_row(dof1);
                auto crd2 = coordinates_row(dof2);
                auto crd3 = coordinates_row(dof3);
                Kokkos::atomic_add(&stress_row(0),
                                   static_cast<Scalar>(crd1 * f1));
                Kokkos::atomic_add(&stress_row(1),
                                   static_cast<Scalar>(crd2 * f2));
                Kokkos::atomic_add(&stress_row(2),
                                   static_cast<Scalar>(crd3 * f3));
                Kokkos::atomic_add(
                    &stress_row(3),
                    static_cast<Scalar>(0.5 * (crd2 * f3 + crd3 * f2)));
                Kokkos::atomic_add(
                    &stress_row(4),
                    static_cast<Scalar>(0.5 * (crd1 * f3 + crd3 * f1)));
                Kokkos::atomic_add(
                    &stress_row(5),
                    static_cast<Scalar>(0.5 * (crd1 * f2 + crd2 * f1)));
              });
          team_member.team_barrier();
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team_member, 6), [=](const int i) {
                stress_row(i) *= scale_factor_d(team_member.league_rank()) /
                                 volume_d(team_member.league_rank());
              });
        });
    stress.template modify<ExeSpace>();
  }
  template <typename ExeSpace, typename T1, typename T2>
  void updateVolume(T1 deformation_gradients, T2 current_volume)
  {
    deformation_gradients.template sync<ExeSpace>();
    current_volume.template sync<ExeSpace>();
    current_volume.template modify<ExeSpace>();
    auto deformation_gradients_d =
        deformation_gradients.template view<ExeSpace>();
    auto current_volume_d = current_volume.template view<ExeSpace>();
    if (deformation_gradients.extent(0) != current_volume.extent(0))
    {
      std::cerr << "The volume and deformation gradients must have the same "
                   "extent corresponding to the number of RVEs."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExeSpace>(0, deformation_gradients_d.extent(0)),
        KOKKOS_LAMBDA(const int i) {
          auto dg = Kokkos::subview(deformation_gradients_d, i, Kokkos::ALL(),
                                    Kokkos::ALL());
          // the jacobian
          auto J = dg(0, 0) * (dg(1, 1) * dg(2, 2) - dg(1, 2) * dg(2, 1)) -
                   dg(0, 1) * (dg(1, 0) * dg(2, 2) - dg(1, 2) * dg(2, 0)) +
                   dg(0, 2) * (dg(1, 0) * dg(2, 1) - dg(1, 1) * dg(2, 0));
          current_volume_d(i) = J * current_volume_d(i);
        });
  }
  template <typename Scalar>
  KOKKOS_FORCEINLINE_FUNCTION Scalar
  getLinearReactionForce(Scalar orig_length,
                         Scalar length,
                         Scalar elastic_modulus,
                         Scalar area
                         )
  {
    // abaqus ...
    // double length_ratio = length / orig_length;
    // double log_strain = log(length_ratio);
    // return elastic_modulus * area * log_strain / length_ratio;
    // Normal 
    Scalar length_ratio = length / orig_length;
    Scalar green_strain = 1.0 / 2.0 * (length_ratio * length_ratio - 1);
    // elastic_modulus * area
    return length_ratio * elastic_modulus * area * green_strain;
  }
  template <typename Derived,
            typename Scalar,
            typename Ordinal,
            typename ExeSpace>
  struct BaseBatchedExplicit
  {
    using PST = PackedData<Scalar*, ExeSpace>;
    using POT = PackedData<Ordinal*, ExeSpace>;
    using ConnectivityType = PackedData<Ordinal*[2], ExeSpace>;
    using HostMemorySpace = typename PST::host_mirror_space;
    using DeviceMemorySpace = typename PST::memory_space;
    // read only view types
    using ROSV = Kokkos::View<const Scalar *, ExeSpace>;
    using ROOV = Kokkos::View<const Ordinal *, ExeSpace>;
    using ConnectivityViewType = Kokkos::View<const Ordinal*[2]>;
    // Random Access View types
    using RASV = ROSV;  // Kokkos::View<const Scalar *,
                        // ExeSpace,
                        // Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
    using RAOV = ROOV;  // Kokkos::View<const Ordinal *,
                        // ExeSpace,
    // Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
    // read write view types
    using RWSV = Kokkos::View<Scalar *, ExeSpace>;
    using RWOV = Kokkos::View<Ordinal *, ExeSpace>;
    static bool run(int num_rves,
                    ConnectivityType connectivity,
                    PST original_coordinates,
                    PST current_coordinates,
                    PST displacement,
                    PST velocity,
                    PST acceleration,
                    PST force_total,
                    PST force_internal,
                    PST force_external,
                    PST force_damping,
                    PST nodal_mass,
                    PST original_length,
                    PST current_length,
                    PST residual,
                    PST fiber_elastic_modulus,
                    PST fiber_area,
                    PST fiber_density,
                    PST viscous_damping_coefficient,
                    PST critical_time_scale_factor,
                    POT displacement_boundary_dof)
    {
      return Derived::run(
          num_rves, connectivity, original_coordinates, current_coordinates,
          displacement, velocity, acceleration, force_total, force_internal,
          force_external, force_damping, nodal_mass, original_length,
          current_length, residual, fiber_elastic_modulus, fiber_area,
          fiber_density, viscous_damping_coefficient,
          critical_time_scale_factor, displacement_boundary_dof);
    }

    protected:
    KOKKOS_FORCEINLINE_FUNCTION
    static void getCurrentCoord(const Ordinal i,
                                ROSV coords,
                                ROSV u,
                                RWSV current_coords)
    {
      current_coords(i) = coords(i) + u(i);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void getElementLength(const Ordinal i,
                                 ROSV coords,
                                 ConnectivityViewType connectivity,
                                 RWSV l0)
    {
      Ordinal n1 = connectivity(i,0)*3;
      Ordinal n2 = connectivity(i,1)*3;
      Scalar x1 = coords(n2) - coords(n1);
      Scalar x2 = coords(n2 + 1) - coords(n1 + 1);
      Scalar x3 = coords(n2 + 2) - coords(n1 + 2);
      l0(i) = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void getForceLoop1(const Ordinal i,
                              RWSV f_int)
    {
      f_int(i) = 0;
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void getForceLoop2(const Ordinal i,
                                ROSV l0,
                                ROSV l,
                                RASV current_coords,
                                ConnectivityViewType connectivity,
                                RWSV f_int,
                                Scalar elastic_modulus,
                                Scalar area)
    {
      Scalar local_l = l(i);
      Scalar frc = getLinearReactionForce(l0(i), local_l, elastic_modulus,area);
      Scalar frc_ovr_l = frc/local_l;
      auto n1 = connectivity(i,0)*3;
      auto n2 = connectivity(i,1)*3;
      Scalar elem_nrm_1 =
          (current_coords(n2) - current_coords(n1)) *frc_ovr_l;
      Scalar elem_nrm_2 =
          (current_coords(n2 + 1) - current_coords(n1 + 1)) *frc_ovr_l;
      Scalar elem_nrm_3 =
          (current_coords(n2 + 2) - current_coords(n1 + 2)) *frc_ovr_l;
      Kokkos::atomic_sub(&f_int(n1), elem_nrm_1);
      Kokkos::atomic_add(&f_int(n2), elem_nrm_1);
      Kokkos::atomic_sub(&f_int(n1 + 1), elem_nrm_2);
      Kokkos::atomic_add(&f_int(n2 + 1), elem_nrm_2);
      Kokkos::atomic_sub(&f_int(n1 + 2), elem_nrm_3);
      Kokkos::atomic_add(&f_int(n2 + 2), elem_nrm_3);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static Scalar getForceLoop3(const Ordinal i,
                                ROSV mass_matrix,
                                ROSV v,
                                ROSV f_int,
                                RWSV f,
                                Scalar viscous_damping_coefficient)
    {
      // viscous damping*constant*mass matrix * velocity
      Scalar f_damp = viscous_damping_coefficient * mass_matrix(i / 3) * v(i);
      // here we assume that the external force is zero on all free
      // degrees of freedom
      Scalar local_residual = -f_int(i);
      f(i) = local_residual - f_damp;
      return local_residual * local_residual;
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void updateAcceleration(const Ordinal i,
                                   ROSV mass_matrix,
                                   ROSV f,
                                   RWSV a)
    {
      a(i) = (1.0 / mass_matrix(i / 3)) * f(i);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void update(const Ordinal i, ROSV a, Scalar dt, RWSV v)
    {
      Kokkos::atomic_add(&v(i), dt * a(i));
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void updateAccelVel(const Ordinal i,
                               Scalar dt,
                               ROSV mass_matrix,
                               ROSV f,
                               RWSV v)
    {
      Scalar a_local = (1.0 / mass_matrix(i / 3)) * f(i);
      Kokkos::atomic_add(&v(i), dt * a_local);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void applyBC(const Ordinal i,
                        Scalar amp_t,
                        ROOV dof,
                        ROSV values,
                        RWSV u)
    {
      u(dof(i)) = values(i) * amp_t;
    }
  };
}  // namespace bio
#endif
