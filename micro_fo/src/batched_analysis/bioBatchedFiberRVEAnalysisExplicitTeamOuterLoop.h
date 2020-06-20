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
        //TeamOuterLoop::getElementLength(i, coords, connectivity, l0);
      Ordinal n1 = connectivity(i * 2)*3;
      Ordinal n2 = connectivity(i * 2 + 1)*3;
      Scalar x1 = coords(n2) - coords(n1);
      Scalar x2 = coords(n2 + 1) - coords(n1 + 1);
      Scalar x3 = coords(n2 + 2) - coords(n1 + 2);
      Scalar tmp =x1 * x1 + x2 * x2 + x3 * x3;
      l0(i) = sqrt(tmp);
      });
    }
    template <typename ExePolicy, typename T, typename T2>
    KOKKOS_INLINE_FUNCTION static void getForces(T team_member,
                                                   ExePolicy element_policy,
                                                   ExePolicy dof_policy,
                                                   T2 scratch,
                                                   ROSV l0,
                                                   ROSV l,
                                                   ROSV v,
                                                   RASV current_coords,
                                                   ROOV connectivity,
                                                   ROSV mass_matrix,
                                                   RWSV f_int,
                                                   RWSV f)
    {
      // swap force arrays and zero internal forces
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::getForceLoop1(i, f_int);
      });
      team_member.team_barrier();
      Kokkos::parallel_for(
          element_policy,
          [=](const Ordinal i) {
            TeamOuterLoop::getForceLoop2(
                i, scratch, l0, l,
                current_coords, connectivity,  f_int);
          });
      // --expt-relaxed-constexpr is needed for this functionality
      // to work properly. Here we find the minimum element size
      // so we can find the overall dt
      Scalar dt = std::numeric_limits<Scalar>::max();
      Kokkos::Min<Scalar> min_reducer(dt);
      Kokkos::parallel_reduce(
          element_policy,
          [=](const Ordinal i, Scalar & dt_crit_elem) {
            min_reducer.join(dt_crit_elem, l(i));
          },
          min_reducer);
      Kokkos::single(Kokkos::PerTeam(team_member), [=]()
      {
      // sqrt(elastic modulus/fiber density
        Scalar sound_speed = sqrt(scratch(1) / scratch(3));
        scratch(4) = dt/sound_speed;
      });
      // barrier needed between loop 2 and loop3
      team_member.team_barrier();

      Scalar residual = 0;
      Kokkos::parallel_reduce(
          dof_policy,
          [=](const Ordinal i, Scalar & residual_update) {
            residual_update +=
                TeamOuterLoop::getForceLoop3(i, scratch, mass_matrix, v, f_int, f);
          },
          residual);
      Kokkos::single(Kokkos::PerTeam(team_member), [=]()
      {
        scratch(5) = sqrt(residual);
      });
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
                                                       RWSV v)
    {
      // 2*ndof loads, 2*ndof writes
      // 2*ndof multiply, 2*ndof divide
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
        TeamOuterLoop::updateAccelVel(i, dt, mass_matrix, f, v);
      });
    }
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void updates(ExePolicy dof_policy,
                                               ROSV velocity,
                                               ROSV mass_matrix,
                                               ROSV f,
                                               Scalar dt,
                                               RWSV current_coordinates)
                                               //RWSV displacement)
    {
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
         
         Scalar a_local = (1/ mass_matrix(i / 3)) * f(i);
         Scalar v_local = velocity(i) + dt * a_local;
         current_coordinates(i) = current_coordinates(i) + dt*v_local;
         //displacement(i) = dt*v_local+displacement(i);
         //Kokkos::atomic_add(&v(i), dt * a_local);
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
      
      /*
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
      */

    using ScratchPad =
        Kokkos::View<Scalar*,
                     typename ExeSpace::scratch_memory_space,
                     Kokkos::MemoryUnmanaged>;
      constexpr int num_scratch = 10;
      auto shared_bytes = ScratchPad::shmem_size(num_scratch);

      // run the main loop
      Kokkos::parallel_for(
          outer_policy.set_scratch_size(0,Kokkos::PerTeam(shared_bytes)),
          KOKKOS_LAMBDA(member_type team_member) {
            ScratchPad scratch(team_member.team_scratch(0), num_rves*num_scratch);
            Ordinal i = team_member.league_rank();
            //// get the device views of the data
            auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
            auto current_coordinates_row =
                current_coordinates.template getRow<ExeSpace>(i);
            auto displacement_row = displacement.template getRow<ExeSpace>(i);
            auto velocity_row = velocity.template getRow<ExeSpace>(i);
            auto force_total_row = force_total.template getRow<ExeSpace>(i);
            auto force_internal_row =
                force_internal.template getRow<ExeSpace>(i);
            auto nodal_mass_row = nodal_mass.template getRow<ExeSpace>(i);
            auto original_length_row =
                original_length.template getRow<ExeSpace>(i);
            auto current_length_row =
                current_length.template getRow<ExeSpace>(i);
            auto displacement_boundary_dof_row =
                displacement_boundary_dof.template getRow<ExeSpace>(i);
            // get the scalar values
            //Scalar visc_damp_coeff =
            //    viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
            constexpr Scalar crit_time_scale_factor = 0.8;
            Kokkos::single(Kokkos::PerTeam(team_member), [=]()
            {
                  // viscous damping coefficient
                  scratch(0) = 
                                viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
                  // fiber elastic modulus
                  scratch(1) = 
                                  fiber_elastic_modulus.template getRow<ExeSpace>(i)(0);
                  // fiber_area
                  scratch(2) = 
                                  fiber_area.template getRow<ExeSpace>(i)(0);
                  // fiber_density
                  scratch(3) = 
                                  fiber_density.template getRow<ExeSpace>(i)(0);
                  // critical time step
                  scratch(4) = dt_crit_d(i);
                  // residual
                  scratch(5) = 10.0;
                  // scratch(6) dt_nphalf
              });
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
            //  need ot have a barrier
            //  to make sure any scratch data has come in properly
            team_member.team_barrier();
            do
            {
              // Kokkos::single was causing an extra 5 registers here! (sm_70)
              if(team_member.team_rank() == 0)
              {
                Scalar dt_nphalf = crit_time_scale_factor * scratch(4)*0.5;
                scratch(6) = dt_nphalf;
              }
              team_member.team_barrier();
              TeamOuterLoop::updates(free_dof_policy, velocity_row, nodal_mass_row, force_total_row,
                                     scratch(6), current_coordinates_row);
              team_member.team_barrier();
              TeamOuterLoop::getElementLengths(
                  element_policy, current_coordinates_row, connectivity_row,
                  current_length_row);
              team_member.team_barrier();
              TeamOuterLoop::getForces(
                  team_member, element_policy, free_dof_policy, scratch,
                   original_length_row, current_length_row,
                  velocity_row, current_coordinates_row, connectivity_row,
                  nodal_mass_row, force_internal_row,
                  force_total_row);
              team_member.team_barrier();
              TeamOuterLoop::updateAccelVels(
                  free_dof_policy, scratch(6), nodal_mass_row,
                  force_total_row, velocity_row);
              team_member.team_barrier();
            } while (scratch(5) > 1E-6);
          });
      /*
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
            // Need to fill the displacement row with current_coords-initial_coords
            //
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
      */

      return true;
    }
  };
}  // namespace bio
#endif
