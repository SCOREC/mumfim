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
    template <typename ExePolicy, typename T>
    KOKKOS_INLINE_FUNCTION static void getElementLengths(
        ExePolicy element_policy,
        ROSV coords,
        ROOV connectivity,
        RWSV l0, T scratch)
    {
      Scalar l_min;
      Kokkos::Min<Scalar> min_reducer(l_min);
      Kokkos::parallel_reduce(element_policy, [=](Ordinal i, Scalar& local) {
        TeamOuterLoop::getElementLength(i, coords, connectivity, l0);
        min_reducer.join(local, l0(i));
      }, min_reducer);
      // might need a barrier here...
      if(element_policy.member.team_rank() == 0)
      {
        Scalar sound_speed = sqrt(scratch(1) / scratch(3));
        scratch(4) = l_min/sound_speed;
      }
    }
    template <typename ExePolicy, typename T>
    KOKKOS_INLINE_FUNCTION static void getForces(  ExePolicy element_policy,
                                                   ExePolicy dof_policy,
                                                   T scratch,
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
      // barrier needed between loop 1 and loop 2
      dof_policy.member.team_barrier();
      Kokkos::parallel_for(
          element_policy,
          [=](const Ordinal i) {
            TeamOuterLoop::getForceLoop2(
                i, scratch, l0, l,
                current_coords, connectivity,  f_int);
          });
      // barrier needed between loop 2 and loop3
      dof_policy.member.team_barrier();

      Scalar residual = 0;
      Kokkos::parallel_reduce(
          dof_policy,
          [=](const Ordinal i, Scalar & residual_update) {
            residual_update +=
                TeamOuterLoop::getForceLoop3(i, scratch, mass_matrix, v, f_int, f);
          },
          residual);
      Kokkos::single(Kokkos::PerTeam(dof_policy.member), [=]()
      {
        scratch(5) = sqrt(residual);
      });
      // we don't put barrier here, because we expect the user to place abarrier after the call if needed
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
                                               RWSV velocity,
                                               ROSV mass_matrix,
                                               ROSV f,
                                               Scalar dt,
                                               RWSV current_coordinates)
    {
      Kokkos::parallel_for(dof_policy, [=](const Ordinal i) {
         
         Scalar a_local = (1/ mass_matrix(i / 3)) * f(i);
         Scalar v_local = velocity(i) + dt * a_local;
         velocity(i) = v_local;
         current_coordinates(i) = current_coordinates(i) + 2*dt*v_local;
      });
    }
    // note ever function we call has been hand tuned to reduce regsiter count so we
    // lie at the 64 reg boundary (so we can acheive 50% occupancy). Please verify
    // that you are not incraseing the reg count using the --resource-usage flag
    // when modifying this code!
    static bool run(Ordinal num_rves,
                    POT connectivity,
                    PST original_coordinates,
                    PST current_coordinates,
                    PST displacement,
                    PST velocity,
                    PST force_total,
                    PST force_internal,
                    PST nodal_mass,
                    PST original_length,
                    PST current_length,
                    PST fiber_elastic_modulus,
                    PST fiber_area,
                    PST fiber_density,
                    PST viscous_damping_coefficient,
                    PST critical_time_scale_factor,
                    POT displacement_boundary_dof,
                    Ordinal num_threads)
    {
      using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace, Kokkos::IndexType<Ordinal>, Kokkos::Schedule<Kokkos::Static>>;
      using member_type = typename OuterPolicyType::member_type;
      OuterPolicyType outer_policy(num_rves, num_threads);
      Kokkos::View<Scalar*,ExeSpace> dt_crit_d("dt_crit", num_rves);
      Kokkos::View<Scalar*,ExeSpace> elastic_modulus_d("elastic_modulus", num_rves);
      using ScratchPad =
          Kokkos::View<Scalar*,
                       typename ExeSpace::scratch_memory_space,
                       Kokkos::MemoryUnmanaged>;
        constexpr int num_scratch = 10;
        auto shared_bytes = ScratchPad::shmem_size(num_scratch);
      
      // initialize data
      // This finds the first dt crit, and initializes the current coordinates
      // to the original coordinates+u
      Kokkos::parallel_for("RunInitialization",
          outer_policy.set_scratch_size(0,Kokkos::PerTeam(shared_bytes)),
          KOKKOS_LAMBDA(member_type team_member) {
            ScratchPad scratch(team_member.team_scratch(0), num_scratch);
            Ordinal i = team_member.league_rank();
            //// get the device views of the data
            auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
            auto original_coordinates_row =
                original_coordinates.template getRow<ExeSpace>(i);
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
            //Kokkos::single(Kokkos::PerTeam(team_member), [=]()
            // Using kokkos single uses ~3 registers here (sm_70)
            if(team_member.team_rank() == 0)
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
                  scratch(4) = 0;
                  // residual
                  scratch(5) = 10.0;
                  // scratch(6) dt_nphalf
                  // scratch(7) l_min
              }
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
            // need to barrier here since the remaining policies need to use the scratch
            // data
            team_member.team_barrier();
            TeamOuterLoop::getElementLengths(
                element_policy, original_coordinates_row, connectivity_row,
                original_length_row, scratch);
            // very important that we update the coordinates with the initial displacements
            // the main loop only updates coordinates by using c(n+1) = c(n) + du
            TeamOuterLoop::getCurrentCoords(
                dof_policy, original_coordinates_row, displacement_row,
                current_coordinates_row);
            team_member.team_barrier();
            TeamOuterLoop::getElementLengths(
                element_policy, current_coordinates_row, connectivity_row,
                current_length_row, scratch);
            team_member.team_barrier();
            // presumably, we only need to do this over the free DOFs
            TeamOuterLoop::getForces(
                element_policy, dof_policy, scratch,
                original_length_row, current_length_row,
                velocity_row, current_coordinates_row, connectivity_row,
                nodal_mass_row, force_internal_row,
                force_total_row);
            team_member.team_barrier();
             dt_crit_d(i) = scratch(4);
          });
      Kokkos::fence();
      // run the main loop
      Kokkos::parallel_for(
          outer_policy.set_scratch_size(0,Kokkos::PerTeam(shared_bytes)),
          KOKKOS_LAMBDA(member_type team_member) {
            ScratchPad scratch(team_member.team_scratch(0), num_scratch);
            Ordinal i = team_member.league_rank();
            //// get the device views of the data
            auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
            auto current_coordinates_row =
                current_coordinates.template getRow<ExeSpace>(i);
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
            // we never actually change this, so save some regs and make it constexpr
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
                team_member, 0, force_total_row.extent(0));
            auto free_dof_policy = Kokkos::TeamThreadRange(
                team_member, 0,
                force_total_row.extent(0) -
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
              //if(team_member.team_rank() == 0)
                //printf("Start While %e, %e\n", scratch(4), scratch(5));
              TeamOuterLoop::updates(free_dof_policy, velocity_row, nodal_mass_row, force_total_row,
                                     scratch(6), current_coordinates_row);
              team_member.team_barrier();
              // this policy also computes the critical dt (min element length/speed of sound
              TeamOuterLoop::getElementLengths(
                  element_policy, current_coordinates_row, connectivity_row,
                  current_length_row, scratch);
              team_member.team_barrier();
              TeamOuterLoop::getForces(
                   element_policy, free_dof_policy, scratch,
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
      Kokkos::fence();
      // Finalize Data
      Kokkos::parallel_for("FinalizeData",
          outer_policy.set_scratch_size(0,Kokkos::PerTeam(shared_bytes)),
          KOKKOS_LAMBDA(member_type team_member) {
            ScratchPad scratch(team_member.team_scratch(0), num_scratch);
            Ordinal i = team_member.league_rank();
            //// get the device views of the data
            auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
            auto original_coordinates_row =
                original_coordinates.template getRow<ExeSpace>(i);
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
            // barrier needed to make sure we wait for shared data
            team_member.team_barrier();
            // This extra call is required to getForces with the full dof_policy
            // to fill in the internal forces on the boundary (which is needed
            // for computation of the cauchy stresses).
              TeamOuterLoop::getForces(
                   element_policy, dof_policy, scratch,
                   original_length_row, current_length_row,
                  velocity_row, current_coordinates_row, connectivity_row,
                  nodal_mass_row, force_internal_row,
                  force_total_row);
            // Need to fill the displacement row with current_coords-initial_coords
              Kokkos::parallel_for(dof_policy, [=](const int i){
                  displacement_row(i) = current_coordinates_row(i) - original_coordinates_row(i);
              });
            team_member.team_barrier();
          });

      return true;
    }
  };
}  // namespace bio
#endif
