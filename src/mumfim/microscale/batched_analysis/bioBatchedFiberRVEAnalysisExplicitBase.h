#ifndef MUMFIM_BATCHED_FIBER_RVE_ANALYSIS_BASE_H
#define MUMFIM_BATCHED_FIBER_RVE_ANALYSIS_BASE_H
#include <Kokkos_Core.hpp>
#include <array>
#include <limits>
#include "bioPackedData.h"
#include "bioRVE.h"
//#define TEAM_SIZE 512
#define TEAM_SIZE Kokkos::AUTO
namespace mumfim
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
                                      displacements_row.extent(0)),
              [&](int i) {
                for (int k = 0; k < 3; ++k)
                {
                  displacements_row(i, k) =
                      (deformation_gradient_row(k, 0) - (k == 0 ? 1 : 0)) *
                          coordinates_row(i,0) +
                      (deformation_gradient_row(k, 1) - (k == 1 ? 1 : 0)) *
                          coordinates_row(i,1) +
                      (deformation_gradient_row(k, 2) - (k == 2 ? 1 : 0)) *
                          coordinates_row(i,2) +
                      displacements_row(i, k);
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
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team_member, boundary_dof_row.extent(0)),
              [=](const int i) {
                auto vert = boundary_dof_row(i);
                auto f1 = force_row(vert,0);
                auto f2 = force_row(vert,1);
                auto f3 = force_row(vert,2);
                auto crd1 = coordinates_row(vert,0);
                auto crd2 = coordinates_row(vert,1);
                auto crd3 = coordinates_row(vert,2);
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
    using ConnectivityViewType = Kokkos::View<Ordinal*[2], ExeSpace>;
    using ScalarDofViewType = Kokkos::View<Scalar*[3], ExeSpace>;
    using ConstScalarDofViewType = Kokkos::View<const Scalar*[3], ExeSpace>;
    using ScalarElementViewType = Kokkos::View<Scalar*, ExeSpace>;
    using ConstScalarElementViewType = Kokkos::View<const Scalar*, ExeSpace>;
    using ScalarVertViewType = Kokkos::View<Scalar*, ExeSpace>;
    using ConstScalarVertViewType = Kokkos::View<const Scalar*, ExeSpace>;

    protected:
    KOKKOS_FORCEINLINE_FUNCTION
    static void getCurrentCoord(const Ordinal i,
                                ConstScalarDofViewType coords,
                                ConstScalarDofViewType u,
                                ScalarDofViewType current_coords)
    {
      for(int j=0; j<3; ++j)
      {
        current_coords(i,j) = coords(i,j) + u(i,j);
      }
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void getElementLength(const Ordinal i,
                                 ConstScalarDofViewType coords,
                                 ConnectivityViewType connectivity,
                                 ScalarElementViewType l0)
    {
      Ordinal n1 = connectivity(i,0);
      Ordinal n2 = connectivity(i,1);
      Scalar x1 = coords(n2,0) - coords(n1,0);
      Scalar x2 = coords(n2,1) - coords(n1,1);
      Scalar x3 = coords(n2,2) - coords(n1,2);
      l0(i) = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void getForceLoop1(const Ordinal i,
                              ScalarDofViewType f_int)
    {
      for(int j=0; j<3; ++j)
      {
        f_int(i,j) = 0;
      }
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void getForceLoop2(const Ordinal i,
                                ConstScalarElementViewType l0,
                                ConstScalarElementViewType l,
                                ConstScalarDofViewType current_coords,
                                ConnectivityViewType connectivity,
                                ScalarDofViewType f_int,
                                Scalar elastic_modulus,
                                Scalar area)
    {
      Scalar local_l = l(i);
      Scalar frc = getLinearReactionForce(l0(i), local_l, elastic_modulus,area);
      Scalar frc_ovr_l = frc/local_l;
      auto n1 = connectivity(i,0);
      auto n2 = connectivity(i,1);
      Scalar elem_nrm_1 =
          (current_coords(n2,0) - current_coords(n1,0)) *frc_ovr_l;
      Scalar elem_nrm_2 =
          (current_coords(n2,1) - current_coords(n1,1)) *frc_ovr_l;
      Scalar elem_nrm_3 =
          (current_coords(n2,2) - current_coords(n1,2)) *frc_ovr_l;
      Kokkos::atomic_sub(&f_int(n1,0), elem_nrm_1);
      Kokkos::atomic_add(&f_int(n2,0), elem_nrm_1);
      Kokkos::atomic_sub(&f_int(n1,1), elem_nrm_2);
      Kokkos::atomic_add(&f_int(n2,1), elem_nrm_2);
      Kokkos::atomic_sub(&f_int(n1,2), elem_nrm_3);
      Kokkos::atomic_add(&f_int(n2,2), elem_nrm_3);
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static Scalar getForceLoop3(const Ordinal i,
                                ConstScalarVertViewType mass_matrix,
                                ConstScalarDofViewType v,
                                ConstScalarDofViewType f_int,
                                ScalarDofViewType f,
                                Scalar viscous_damping_coefficient)
    {
      Scalar mass = mass_matrix(i)*viscous_damping_coefficient;
      Scalar residual = 0;
      for(int j=0; j<3; ++j)
      {
        // viscous damping*constant*mass matrix * velocity
        Scalar f_damp = mass*v(i,j);
        // here we assume that the external force is zero on all free
        // degrees of freedom
        Scalar local_residual = -f_int(i,j);
        f(i,j) = local_residual - f_damp;
        residual += local_residual*local_residual;
      }
      return residual;
    }
    KOKKOS_FORCEINLINE_FUNCTION
    static void updateAccelVel(const Ordinal i,
                               Scalar dt,
                               ConstScalarVertViewType mass_matrix,
                               ConstScalarDofViewType f,
                               ScalarDofViewType v)
    {
      Scalar mass = mass_matrix(i);
      for(int j=0; j<3; ++j)
      {
        v(i,j) += dt * f(i,j)/mass;
      }
    }
  };
}  // namespace mumfim
#endif
