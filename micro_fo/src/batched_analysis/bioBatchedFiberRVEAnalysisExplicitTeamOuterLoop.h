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
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static Scalar getForces(ExePolicy element_policy,
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
      residual = 0;
      // swap force arrays and zero internal forces
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::getForceLoop1(i, visc_damp_coeff, mass_matrix, v, f_int,
                                     f_damp);
      });
      // set the internal forces
      Kokkos::Min<Scalar> min_reducer(dt);
      Kokkos::parallel_reduce(
          element_policy,
          [=](const Ordinal i, Scalar & dt_crit_elem) {
            auto local_l = TeamOuterLoop::getForceLoop2(
                i, fiber_elastic_modulus, fiber_area, fiber_density, l0, l,
                current_coords, connectivity, visc_damp_coeff, f_int);
            min_reducer.join(dt_crit_elem, local_l / sound_speed);
          },
          min_reducer);
      residual = 0;
      Kokkos::parallel_reduce(
          dof_policy,
          [=](const Ordinal i, Scalar & residual_update) {
            residual_update +=
                TeamOuterLoop::getForceLoop3(i, f_ext, f_int, f_damp, f);
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
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION
        // FIXME Add init values back!
        static void
        applyBCs(ExePolicy fixed_dof_policy,
                 Scalar amp_t,
                 ROOV dof,
                 ROSV values,
                 RWSV u)
    {
      Kokkos::parallel_for(fixed_dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::applyBC(i, amp_t, dof, values, u);
      });
    }
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void applyAccelVelBCs(
        ExePolicy fixed_dof_policy,
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
      // 4*nfixed reads, 5*nfixed writes
      // 5*nfixed multiply, 2*nfixed adds, 1 divide
      Kokkos::parallel_for(fixed_dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::applyAccelVelBC(i, a_amp, v_amp, visc_damp_coeff, dof,
                                       values, mass_matrix, a, v, f_int, f_ext,
                                       f_damp, f);
      });
    }
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
      using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace>;
      using member_type = typename OuterPolicyType::member_type;
      // OuterPolicyType outer_policy(rves.size(),
      //                             Kokkos::AUTO());
      OuterPolicyType outer_policy(rves.size(), TEAM_SIZE);
      // OuterPolicyType outer_policy(rves.size(), 256);
      Kokkos::parallel_for(
          outer_policy, KOKKOS_LAMBDA(member_type team_member) {
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
            auto residual_row = residual.template getRow<ExeSpace>(i);
            auto displacement_boundary_dof_row =
                displacement_boundary_dof.template getRow<ExeSpace>(i);
            auto displacement_boundary_values_row =
                displacement_boundary_values.template getRow<ExeSpace>(i);

            //Kokkos::View<const Scalar*, ExeSpace, Kokkos::MemoryRandomAccess> test = 
            //    displacement_boundary_values.template getRow<ExeSpace>(i);
            //printf("Random access array %f", test(0));
            // get the scalar values
            Scalar visc_damp_coeff =
                viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
            Scalar crit_time_scale_factor =
                critical_time_scale_factor.template getRow<ExeSpace>(i)(0);
            Scalar elastic_modulus =
                fiber_elastic_modulus.template getRow<ExeSpace>(i)(0);
            Scalar area = fiber_area.template getRow<ExeSpace>(i)(0);
            Scalar density = fiber_density.template getRow<ExeSpace>(i)(0);
            // create the loop policies based on the lengths of the various
            // arrays
            auto element_policy = Kokkos::TeamThreadRange(
                team_member, current_length_row.extent(0));
            auto dof_policy = Kokkos::TeamThreadRange(
                team_member, displacement_row.extent(0));
            auto fixed_dof_policy = Kokkos::TeamThreadRange(
                team_member, displacement_boundary_dof_row.extent(0));
            Scalar residual = 10.0;
            Scalar total_time = 10.0;
            Scalar current_time = 10.0;
            long int n_step = 0;
            Scalar dt_nphalf, t_npone, t_nphalf;
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
            auto dt_crit = TeamOuterLoop::getForces(
                element_policy, dof_policy, elastic_modulus, area, density,
                original_length_row, current_length_row, velocity_row,
                current_coordinates_row, connectivity_row, nodal_mass_row,
                visc_damp_coeff, force_internal_row, force_external_row,
                force_damping_row, force_total_row, residual);
            team_member.team_barrier();
            TeamOuterLoop::updateAccelerations(
                dof_policy, nodal_mass_row, force_total_row, acceleration_row);
            team_member.team_barrier();
            do
            {
              dt_nphalf = crit_time_scale_factor * dt_crit;
              t_npone = current_time + dt_nphalf;
              t_nphalf = 0.5 * (current_time + t_npone);
              assert(dt_nphalf != t_npone != t_nphalf);
              TeamOuterLoop::updates(dof_policy, acceleration_row,
                                     (t_nphalf - current_time), velocity_row);
              team_member.team_barrier();
              TeamOuterLoop::updates(dof_policy, velocity_row, dt_nphalf,
                                     displacement_row);
              team_member.team_barrier();
              TeamOuterLoop::applyBCs(
                  fixed_dof_policy, 1.0, displacement_boundary_dof_row,
                  displacement_boundary_values_row, displacement_row);
              team_member.team_barrier();
              TeamOuterLoop::getCurrentCoords(
                  dof_policy, original_coordinates_row, displacement_row,
                  current_coordinates_row);
              team_member.team_barrier();
              TeamOuterLoop::getElementLengths(
                  element_policy, current_coordinates_row, connectivity_row,
                  current_length_row);
              team_member.team_barrier();
              dt_crit = TeamOuterLoop::getForces(
                  element_policy, dof_policy, elastic_modulus, area, density,
                  original_length_row, current_length_row, velocity_row,
                  current_coordinates_row, connectivity_row, nodal_mass_row,
                  visc_damp_coeff, force_internal_row, force_external_row,
                  force_damping_row, force_total_row, residual);
              team_member.team_barrier();
              TeamOuterLoop::updateAccelVels(dof_policy, t_npone - t_nphalf,
                                             nodal_mass_row, force_total_row,
                                             acceleration_row, velocity_row);
              team_member.team_barrier();
              // the derivative and second derivative amplitudes are zero
              // since we apply the displacement boundary condition all at once
              TeamOuterLoop::applyAccelVelBCs(
                  fixed_dof_policy, (Scalar)0, (Scalar)0, visc_damp_coeff,
                  displacement_boundary_dof_row,
                  displacement_boundary_values_row, nodal_mass_row,
                  acceleration_row, velocity_row, force_internal_row,
                  force_external_row, force_damping_row, force_total_row);
              current_time = t_npone;
              ++n_step;
              team_member.team_barrier();
              // if(team_member.league_rank() == 0 && team_member.team_rank() ==
              // 0)
              //  printf("%ld, %f\n", n_step, residual);
              //} while ((current_time < total_time) || (residual > 1E-6));
            } while ((current_time < total_time) || (residual > 1E-4));
            if (team_member.league_rank() == 0 && team_member.team_rank() == 0)
              printf("%ld, %f\n", n_step, residual);
          });
      return true;
    }
  };
}  // namespace bio
#endif
