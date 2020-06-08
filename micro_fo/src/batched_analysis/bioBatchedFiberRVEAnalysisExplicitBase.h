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
            typename T,
            typename T2,
            typename T3,
            int DIM = 3>
  void applyIncrementalDeformationToDisplacement(T deformation_gradient,
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
    OuterPolicyType outer_policy(deformation_gradient.extent(0), TEAM_SIZE);
    Kokkos::parallel_for(
        "applyIncrementalDfm", outer_policy,
        KOKKOS_LAMBDA(member_type team_member) {
          int rve = team_member.league_rank();
          auto displacements_row = displacements.template getRow<ExeSpace>(rve);
          auto coordinates_row = coordinates.template getRow<ExeSpace>(rve);
          auto deformation_gradient_row = Kokkos::subview(
              deformation_gradient_d, rve, Kokkos::ALL(), Kokkos::ALL());
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
            int DIM = 3>
  void applyIncrementalDeformationToBoundary(T deformation_gradient,
                                             T2 coordinates,
                                             T3 boundary_dof,
                                             T2 boundary_values)
  {
    deformation_gradient.template sync<ExeSpace>();
    auto deformation_gradient_d =
        deformation_gradient.template view<ExeSpace>();
    coordinates.template sync<ExeSpace>();
    boundary_dof.template sync<ExeSpace>();
    boundary_values.template sync<ExeSpace>();
    boundary_values.template modify<ExeSpace>();
    using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace>;
    using member_type = typename OuterPolicyType::member_type;
    OuterPolicyType outer_policy(deformation_gradient.extent(0), TEAM_SIZE);
    Kokkos::parallel_for(
        "applyIncrementalDfm", outer_policy,
        KOKKOS_LAMBDA(member_type team_member) {
          int rve = team_member.league_rank();
          auto deformation_gradient_row = Kokkos::subview(
              deformation_gradient_d, rve, Kokkos::ALL(), Kokkos::ALL());
          auto coordinates_row = coordinates.template getRow<ExeSpace>(rve);
          auto boundary_dof_row = boundary_dof.template getRow<ExeSpace>(rve);
          auto boundary_values_row =
              boundary_values.template getRow<ExeSpace>(rve);
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team_member,
                                      boundary_dof_row.extent(0) / DIM),
              [&](int i) {
                auto dof1 = boundary_dof_row(i * DIM);
                auto dof2 = boundary_dof_row(i * DIM + 1);
                auto dof3 = boundary_dof_row(i * DIM + 2);
                for (int k = 0; k < 3; ++k)
                {
                  boundary_values_row(i * DIM + k) =
                      (deformation_gradient_row(k, 0) - (k == 0 ? 1 : 0)) *
                          coordinates_row(dof1) +
                      (deformation_gradient_row(k, 1) - (k == 1 ? 1 : 0)) *
                          coordinates_row(dof2) +
                      (deformation_gradient_row(k, 2) - (k == 2 ? 1 : 0)) *
                          coordinates_row(dof3) +
                      boundary_values_row(i * DIM + k);
                }
              });
        });
  }
  template <typename ExeSpace, typename T, typename T2, typename T3>
  void computeCauchyStress(T2 boundary_dofs, T3 coordinates, T3 force, T stress)
  {
    using Scalar = typename T::value_type;
    boundary_dofs.template sync<ExeSpace>();
    coordinates.template sync<ExeSpace>();
    force.template sync<ExeSpace>();
    stress.template modify<ExeSpace>();
    auto stress_d = stress.template view<ExeSpace>();
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
        });
  }
  template <typename Scalar>
  KOKKOS_FORCEINLINE_FUNCTION Scalar
  getLinearReactionForce(Scalar orig_length,
                         Scalar length,
                         Scalar elastic_modulus,
                         Scalar area)
  {
    // abaqus ...
    // double length_ratio = length / orig_length;
    // double log_strain = log(length_ratio);
    // return elastic_modulus * area * log_strain / length_ratio;
    Scalar length_ratio = length / orig_length;
    Scalar green_strain = 1.0 / 2.0 * (length_ratio * length_ratio - 1);
    return length_ratio * elastic_modulus * area * green_strain;
  }
  template <typename Derived,
            typename Scalar,
            typename Ordinal,
            typename ExeSpace>
  struct BaseBatchedExplicit
  {
    using PST = PackedData<Scalar, Ordinal, ExeSpace>;
    using POT = PackedData<Ordinal, Ordinal, ExeSpace>;
    using HostMemorySpace = typename PST::host_mirror_space;
    using DeviceMemorySpace = typename PST::memory_space;
    // read only view types
    using ROSV = Kokkos::View<const Scalar *, ExeSpace>;
    using ROOV = Kokkos::View<const Ordinal *, ExeSpace>;
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
    static bool run(const std::vector<RVE> & rves,
                    POT connectivity,
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
                    POT displacement_boundary_dof,
                    PST displacement_boundary_values)
    {
      return Derived::run(
          rves, connectivity, original_coordinates, current_coordinates,
          displacement, velocity, acceleration, force_total, force_internal,
          force_external, force_damping, nodal_mass, original_length,
          current_length, residual, fiber_elastic_modulus, fiber_area,
          fiber_density, viscous_damping_coefficient,
          critical_time_scale_factor, displacement_boundary_dof,
          displacement_boundary_values);
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
                                 ROOV connectivity,
                                 RWSV l0)
    {
      Ordinal n1 = connectivity(i * 2);
      Ordinal n2 = connectivity(i * 2 + 1);
      Scalar x1 = coords(n2 * 3) - coords(n1 * 3);
      Scalar x2 = coords(n2 * 3 + 1) - coords(n1 * 3 + 1);
      Scalar x3 = coords(n2 * 3 + 2) - coords(n1 * 3 + 2);
      l0(i) = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void getForceLoop1(const Ordinal i,
                              Scalar visc_damp_coeff,
                              ROSV mass_matrix,
                              ROSV v,
                              RWSV f_int,
                              RWSV f_damp)
    {
      f_int(i) = 0;
      f_damp(i) = visc_damp_coeff * mass_matrix(i / 3) * v(i);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static Scalar getForceLoop2(const Ordinal i,
                                Scalar fiber_elastic_modulus,
                                Scalar fiber_area,
                                Scalar fiber_density,
                                ROSV l0,
                                ROSV l,
                                RASV current_coords,
                                ROOV connectivity,
                                Scalar visc_damp_coeff,
                                RWSV f_int)
    {
      Scalar local_l = l(i);
      // FIXME the force reactions are all messed up! need to figure out
      // how to get the struct data in here?
      Scalar frc = getLinearReactionForce(l0(i), local_l, fiber_elastic_modulus,
                                          fiber_area);
      auto n1 = connectivity(i * 2);
      auto n2 = connectivity(i * 2 + 1);
      Scalar elem_nrm_1 =
          (current_coords(n2 * 3) - current_coords(n1 * 3)) / local_l;
      Scalar elem_nrm_2 =
          (current_coords(n2 * 3 + 1) - current_coords(n1 * 3 + 1)) / local_l;
      Scalar elem_nrm_3 =
          (current_coords(n2 * 3 + 2) - current_coords(n1 * 3 + 2)) / local_l;
      Kokkos::atomic_add(&f_int(n1 * 3),
                         static_cast<Scalar>(-frc * elem_nrm_1));
      Kokkos::atomic_add(&f_int(n1 * 3 + 1),
                         static_cast<Scalar>(-frc * elem_nrm_2));
      Kokkos::atomic_add(&f_int(n1 * 3 + 2),
                         static_cast<Scalar>(-frc * elem_nrm_3));
      Kokkos::atomic_add(&f_int(n2 * 3), static_cast<Scalar>(frc * elem_nrm_1));
      Kokkos::atomic_add(&f_int(n2 * 3 + 1),
                         static_cast<Scalar>(frc * elem_nrm_2));
      Kokkos::atomic_add(&f_int(n2 * 3 + 2),
                         static_cast<Scalar>(frc * elem_nrm_3));
      return local_l;
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static Scalar getForceLoop3(const Ordinal i,
                                ROSV f_ext,
                                ROSV f_int,
                                ROSV f_damp,
                                RWSV f)
    {
      Scalar local_residual = f_ext(i) - f_int(i);
      f(i) = local_residual - f_damp(i);
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
      //v(i) += dt * a(i);
      Kokkos::atomic_add(&v(i), dt * a(i));
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void updateAccelVel(const Ordinal i,
                               Scalar dt,
                               ROSV mass_matrix,
                               ROSV f,
                               RWSV a,
                               RWSV v)
    {
      Scalar a_local = (1.0 / mass_matrix(i / 3)) * f(i);
      a(i) = a_local;
      //v(i) += dt * a_local;
      Kokkos::atomic_add(&v(i), dt*a_local);
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
    KOKKOS_FORCEINLINE_FUNCTION
    static void applyAccelVelBC(const Ordinal i,
                                Scalar a_amp,
                                Scalar v_amp,
                                Scalar visc_damp_coeff,
                                ROOV dof,
                                ROSV values,
                                ROSV mass_matrix,
                                RWSV a,
                                RWSV v,
                                RWSV f_int,
                                RWSV f_ext,
                                RWSV f_damp,
                                RWSV f)
    {
      auto local_dof = dof(i);
      Scalar value = values(i);
      Scalar a_local = value * a_amp;
      Scalar v_local = value * v_amp;
      Scalar mass = mass_matrix(local_dof / 3);
      Scalar f_damp_local = visc_damp_coeff * mass * v_local;
      Scalar f_inertial = mass * a_local;
      a(local_dof) = a_local;
      v(local_dof) = v_local;
      f_damp(local_dof) = f_damp_local;
      f_ext(local_dof) = f_inertial + f_int(local_dof) + f_damp_local;
      f(local_dof) = f_inertial;
    }
  };
}  // namespace bio
#endif
