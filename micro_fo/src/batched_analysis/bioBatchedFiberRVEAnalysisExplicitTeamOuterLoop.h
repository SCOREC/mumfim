#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_TEAM_OUTER_LOOP_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_TEAM_OUTER_LOOP_H__
#include "bioBatchedFiberRVEAnalysisExplicitBase.h"
namespace bio
{
  template <typename Scalar, typename Ordinal, typename ExeSpace>
  struct TeamOuterLoop
      : public BaseBatchedExplicit<TeamOuterLoop<Scalar, Ordinal, ExeSpace>,
                                   Scalar,
                                   Ordinal,
                                   ExeSpace>
  {
    using BaseType =
        BaseBatchedExplicit<TeamOuterLoop<Scalar, Ordinal, ExeSpace>,
                            Scalar,
                            Ordinal,
                            ExeSpace>;
    using typename BaseType::DeviceMemorySpace;
    using typename BaseType::HostMemorySpace;
    using typename BaseType::POT;
    using typename BaseType::PST;
    using typename BaseType::RAOV;
    using typename BaseType::RASV;
    using typename BaseType::ROOV;
    using typename BaseType::ROSV;
    using typename BaseType::RWOV;
    using typename BaseType::RWSV;
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void getCurrentCoords(ExePolicy dof_policy,
                                                        ROSV coords,
                                                        ROSV u,
                                                        RWSV current_coords)
    {
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::getCurrentCoord(i, coords, u, current_coords);
      });
    }
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void getElementLengths(
        ExePolicy element_policy,
        ROSV coords,
        ROOV connectivity,
        RWSV l0)
    {
      Kokkos::parallel_for(element_policy, [=](Ordinal i) {
        TeamOuterLoop::getElementLength(i, coords, connectivity, l0);
      });
    }
    template <typename ExePolicy, typename T>
    KOKKOS_INLINE_FUNCTION static Scalar getForces(T team_member,
                                                   ExePolicy element_policy,
                                                   ExePolicy dof_policy,
                                                   Scalar fiber_elastic_modulus,
                                                   Scalar fiber_area,
                                                   Scalar fiber_density,
                                                   ROSV l0,
                                                   ROSV l,
                                                   ROSV v,
                                                   RASV current_coords,
                                                   ROOV connectivity,
                                                   ROSV mass_matrix,
                                                   Scalar visc_damp_coeff,
                                                   RWSV f_int,
                                                   RWSV f_ext,
                                                   RWSV f_damp,
                                                   RWSV f,
                                                   Scalar & residual)
    {
      // a large number
      Scalar dt = 1000000000;  // std::numeric_limits<Scalar>::max();
      // element normal vectors
      Scalar sound_speed = sqrt(fiber_elastic_modulus / fiber_density);
      // swap force arrays and zero internal forces
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::getForceLoop1(i, f_int);
      });
      team_member.team_barrier();
      // set the internal forces
      Kokkos::Min<Scalar> min_reducer(dt);
      Kokkos::parallel_reduce(
          element_policy,
          [=](const Ordinal i, Scalar & dt_crit_elem) {
            auto local_l = TeamOuterLoop::getForceLoop2(
                i, fiber_elastic_modulus, fiber_area, fiber_density, l0, l,
                current_coords, connectivity,  f_int);
            min_reducer.join(dt_crit_elem, local_l / sound_speed);
          },
          min_reducer);
      team_member.team_barrier();
      Kokkos::parallel_reduce(
          dof_policy,
          [=](const Ordinal i, Scalar & residual_update) {
            residual_update +=
                TeamOuterLoop::getForceLoop3(i, visc_damp_coeff, mass_matrix, v, f_ext, f_int, f);
          },
          residual);
      residual = sqrt(residual);
      return dt;
    }
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void updateAccelerations(ExePolicy dof_policy,
                                                           ROSV mass_matrix,
                                                           ROSV f,
                                                           RWSV a)
    {
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::updateAcceleration(i, mass_matrix, f, a);
      });
    }
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void updateAccelVels(ExePolicy dof_policy,
                                                       Scalar dt,
                                                       ROSV mass_matrix,
                                                       ROSV f,
                                                       RWSV a,
                                                       RWSV v)
    {
      // 2*ndof loads, 2*ndof writes
      // 2*ndof multiply, 2*ndof divide
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::updateAccelVel(i, dt, mass_matrix, f, a, v);
      });
    }
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void updates(ExePolicy dof_policy,
                                               ROSV a,
                                               Scalar dt,
                                               RWSV v)
    {
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::update(i, a, dt, v);
      });
    }
    static bool run(Ordinal num_rves,
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
                    Ordinal num_threads)
    {
      using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace>;
      using member_type = typename OuterPolicyType::member_type;
      OuterPolicyType outer_policy(num_rves, num_threads);
      Kokkos::View<Scalar*,ExeSpace> dt_crit_d("dt_crit", num_rves);
      Kokkos::View<Scalar*,ExeSpace> elastic_modulus_d("elastic_modulus", num_rves);
      
      // initialize data
      Kokkos::parallel_for("RunInitialization",
          outer_policy,
          KOKKOS_LAMBDA(member_type team_member) {
            Ordinal i = team_member.league_rank();
            //// get the device views of the data
            auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
            auto original_coordinates_row =
                original_coordinates.template getRow<ExeSpace>(i);
            auto current_coordinates_row =
                current_coordinates.template getRow<ExeSpace>(i);
            auto displacement_row = displacement.template getRow<ExeSpace>(i);
            auto velocity_row = velocity.template getRow<ExeSpace>(i);
            auto acceleration_row = acceleration.template getRow<ExeSpace>(i);
            auto force_total_row = force_total.template getRow<ExeSpace>(i);
            auto force_internal_row =
                force_internal.template getRow<ExeSpace>(i);
            auto force_external_row =
                force_external.template getRow<ExeSpace>(i);
            auto force_damping_row = force_damping.template getRow<ExeSpace>(i);
            auto nodal_mass_row = nodal_mass.template getRow<ExeSpace>(i);
            auto original_length_row =
                original_length.template getRow<ExeSpace>(i);
            auto current_length_row =
                current_length.template getRow<ExeSpace>(i);
            auto displacement_boundary_dof_row =
                displacement_boundary_dof.template getRow<ExeSpace>(i);
            // get the scalar values
            Scalar visc_damp_coeff =
                viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
            //Scalar crit_time_scale_factor =
            //    critical_time_scale_factor.template getRow<ExeSpace>(i)(0);
            Scalar elastic_modulus =
                fiber_elastic_modulus.template getRow<ExeSpace>(i)(0);
            Scalar area = fiber_area.template getRow<ExeSpace>(i)(0);
            Scalar density = fiber_density.template getRow<ExeSpace>(i)(0);
            // create the loop policies based on the lengths of the various
            // arrays
            auto element_policy = Kokkos::TeamThreadRange(
                team_member, 0, current_length_row.extent(0));
            auto dof_policy = Kokkos::TeamThreadRange(
                team_member, 0, displacement_row.extent(0));
            auto free_dof_policy = Kokkos::TeamThreadRange(
                team_member, 0,
                displacement_row.extent(0) -
                    displacement_boundary_dof_row.extent(0));
            Scalar residual = 10.0;
            TeamOuterLoop::getElementLengths(
                element_policy, original_coordinates_row, connectivity_row,
                original_length_row);
            TeamOuterLoop::getCurrentCoords(
                dof_policy, original_coordinates_row, displacement_row,
                current_coordinates_row);
            team_member.team_barrier();
            TeamOuterLoop::getElementLengths(
                element_policy, current_coordinates_row, connectivity_row,
                current_length_row);
            team_member.team_barrier();
            dt_crit_d(i) = TeamOuterLoop::getForces(
                team_member, element_policy, free_dof_policy, elastic_modulus,
                area, density, original_length_row, current_length_row,
                velocity_row, current_coordinates_row, connectivity_row,
                nodal_mass_row, visc_damp_coeff, force_internal_row,
                force_external_row, force_damping_row, force_total_row,
                residual);
            team_member.team_barrier();
            TeamOuterLoop::updateAccelerations(free_dof_policy, nodal_mass_row,
                                               force_total_row,
                                               acceleration_row);
          });

      // run the main loop
      Kokkos::parallel_for(
          outer_policy,
          KOKKOS_LAMBDA(member_type team_member) {
            Ordinal i = team_member.league_rank();
            //// get the device views of the data
            auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
            auto original_coordinates_row =
                original_coordinates.template getRow<ExeSpace>(i);
            auto current_coordinates_row =
                current_coordinates.template getRow<ExeSpace>(i);
            auto displacement_row = displacement.template getRow<ExeSpace>(i);
            auto velocity_row = velocity.template getRow<ExeSpace>(i);
            auto acceleration_row = acceleration.template getRow<ExeSpace>(i);
            auto force_total_row = force_total.template getRow<ExeSpace>(i);
            auto force_internal_row =
                force_internal.template getRow<ExeSpace>(i);
            auto force_external_row =
                force_external.template getRow<ExeSpace>(i);
            auto force_damping_row = force_damping.template getRow<ExeSpace>(i);
            auto nodal_mass_row = nodal_mass.template getRow<ExeSpace>(i);
            auto original_length_row =
                original_length.template getRow<ExeSpace>(i);
            auto current_length_row =
                current_length.template getRow<ExeSpace>(i);
            auto displacement_boundary_dof_row =
                displacement_boundary_dof.template getRow<ExeSpace>(i);
            // get the scalar values
            Scalar visc_damp_coeff =
                viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
            constexpr Scalar crit_time_scale_factor = 0.8;
            //constexpr Scalar elastic_modulus = 5.0;
            //Scalar elastic_modulus = elastic_modulus_d(i);
            //Scalar crit_time_scale_factor =
            //    critical_time_scale_factor.template getRow<ExeSpace>(i)(0);
            Scalar elastic_modulus =
                fiber_elastic_modulus.template getRow<ExeSpace>(i)(0);
            Scalar area = fiber_area.template getRow<ExeSpace>(i)(0);
            Scalar density = fiber_density.template getRow<ExeSpace>(i)(0);
            // create the loop policies based on the lengths of the various
            // arrays
            auto element_policy = Kokkos::TeamThreadRange(
                team_member, 0, current_length_row.extent(0));
            auto dof_policy = Kokkos::TeamThreadRange(
                team_member, 0, displacement_row.extent(0));
            auto free_dof_policy = Kokkos::TeamThreadRange(
                team_member, 0,
                displacement_row.extent(0) -
                    displacement_boundary_dof_row.extent(0));
            Scalar residual = 10.0;
            Scalar dt_crit = dt_crit_d(i);
            do
            {
              Scalar dt_nphalf = crit_time_scale_factor * dt_crit*0.5;
              TeamOuterLoop::updates(free_dof_policy, acceleration_row,
                                     dt_nphalf, velocity_row);
              team_member.team_barrier();
              TeamOuterLoop::updates(free_dof_policy, velocity_row,dt_nphalf ,
                                     displacement_row);
              team_member.team_barrier();
              TeamOuterLoop::getCurrentCoords(
                  free_dof_policy, original_coordinates_row, displacement_row,
                  current_coordinates_row);
              team_member.team_barrier();
              TeamOuterLoop::getElementLengths(
                  element_policy, current_coordinates_row, connectivity_row,
                  current_length_row);
              team_member.team_barrier();
              dt_crit = TeamOuterLoop::getForces(
                  team_member, element_policy, free_dof_policy, elastic_modulus,
                  area, density, original_length_row, current_length_row,
                  velocity_row, current_coordinates_row, connectivity_row,
                  nodal_mass_row, visc_damp_coeff, force_internal_row,
                  force_external_row, force_damping_row, force_total_row,
                  residual);
              team_member.team_barrier();
              TeamOuterLoop::updateAccelVels(
                  free_dof_policy, dt_nphalf, nodal_mass_row,
                  force_total_row, acceleration_row, velocity_row);
              team_member.team_barrier();
            } while (residual > 1E-6);
          });
      // Finalize Data
      Kokkos::parallel_for("FinalizeData",
          outer_policy,
          KOKKOS_LAMBDA(member_type team_member) {
            Ordinal i = team_member.league_rank();
            //// get the device views of the data
            auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
            auto original_coordinates_row =
                original_coordinates.template getRow<ExeSpace>(i);
            auto current_coordinates_row =
                current_coordinates.template getRow<ExeSpace>(i);
            auto displacement_row = displacement.template getRow<ExeSpace>(i);
            auto velocity_row = velocity.template getRow<ExeSpace>(i);
            auto acceleration_row = acceleration.template getRow<ExeSpace>(i);
            auto force_total_row = force_total.template getRow<ExeSpace>(i);
            auto force_internal_row =
                force_internal.template getRow<ExeSpace>(i);
            auto force_external_row =
                force_external.template getRow<ExeSpace>(i);
            auto force_damping_row = force_damping.template getRow<ExeSpace>(i);
            auto nodal_mass_row = nodal_mass.template getRow<ExeSpace>(i);
            auto original_length_row =
                original_length.template getRow<ExeSpace>(i);
            auto current_length_row =
                current_length.template getRow<ExeSpace>(i);
            auto displacement_boundary_dof_row =
                displacement_boundary_dof.template getRow<ExeSpace>(i);
            // get the scalar values
            Scalar visc_damp_coeff =
                viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
            Scalar elastic_modulus =
                fiber_elastic_modulus.template getRow<ExeSpace>(i)(0);
            Scalar area = fiber_area.template getRow<ExeSpace>(i)(0);
            Scalar density = fiber_density.template getRow<ExeSpace>(i)(0);
            // create the loop policies based on the lengths of the various
            // arrays
            auto element_policy = Kokkos::TeamThreadRange(
                team_member, 0, current_length_row.extent(0));
            auto dof_policy = Kokkos::TeamThreadRange(
                team_member, 0, displacement_row.extent(0));
            auto free_dof_policy = Kokkos::TeamThreadRange(
                team_member, 0,
                displacement_row.extent(0) -
                    displacement_boundary_dof_row.extent(0));
            Scalar residual = 10.0;
            // This extra call is required to getForces with the full dof_policy
            // to fill in the internal forces on the boundary (which is needed
            // for computation of the cauchy stresses).
            TeamOuterLoop::getForces(
                team_member, element_policy, dof_policy, elastic_modulus, area,
                density, original_length_row, current_length_row, velocity_row,
                current_coordinates_row, connectivity_row, nodal_mass_row,
                visc_damp_coeff, force_internal_row, force_external_row,
                force_damping_row, force_total_row, residual);
          });

      return true;
    }
  };
}  // namespace bio
#endif
