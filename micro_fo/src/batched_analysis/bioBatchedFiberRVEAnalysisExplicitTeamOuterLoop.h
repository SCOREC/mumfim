#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_TEAM_OUTER_LOOP_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_TEAM_OUTER_LOOP_H__
#include "bioBatchedFiberRVEAnalysisExplicitBase.h"

#define UPDATE_FREQ 100

namespace bio
{
  template<typename Scalar, typename Ordinal>
  struct ScratchRegister
  {
    Scalar viscous_damping_coefficient; // (0)
    Scalar fiber_elastic_modulus;  // (1)
    Scalar fiber_area; // (2)
    Scalar fiber_density; // (3)
    Scalar dt_crit; // (4)
    Scalar residual; // (5)
    Scalar dt_nphalf; // (6)
    Ordinal step;
  };

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

    using typename BaseType::ConnectivityType;
    using typename BaseType::ConnectivityViewType;
    using typename BaseType::ScalarDofViewType;
    using typename BaseType::ConstScalarDofViewType;
    using typename BaseType::ScalarElementViewType;
    using typename BaseType::ConstScalarElementViewType;
    using typename BaseType::ScalarVertViewType;
    using typename BaseType::ConstScalarVertViewType;

    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void getCurrentCoords(ExePolicy vert_policy,
                                                        ConstScalarDofViewType coords,
                                                        ConstScalarDofViewType u,
                                                        ScalarDofViewType current_coords)
    {
      Kokkos::parallel_for(vert_policy, [=](const Ordinal i) {
        TeamOuterLoop::getCurrentCoord(i, coords, u, current_coords);
      });
    }
    template <typename ExePolicy, typename T, typename Member>
    KOKKOS_INLINE_FUNCTION static void getElementLengths(
        ExePolicy element_policy,
        ConstScalarDofViewType coords,
        ConnectivityViewType connectivity,
        ScalarElementViewType l0, T scratch,
        const Member& team)
    {
      if(scratch(0).step % UPDATE_FREQ == 0)
      {
        Scalar l_min;
        Kokkos::Min<Scalar> min_reducer(l_min);
        Kokkos::parallel_reduce(element_policy, [=](Ordinal i, Scalar& local) {
          TeamOuterLoop::getElementLength(i, coords, connectivity, l0);
          min_reducer.join(local, l0(i));
        }, min_reducer);
        // might need a barrier here...
        if(team.team_rank() == 0)
        {
          Scalar sound_speed = sqrt(scratch(0).fiber_elastic_modulus / scratch(0).fiber_density);
          scratch(0).dt_crit = l_min/sound_speed;
        }
      }
      else
      {
        Kokkos::parallel_for(element_policy, [=](Ordinal i) {
          TeamOuterLoop::getElementLength(i, coords, connectivity, l0);
        });
      }
    }
    template <typename ExePolicy, typename T, typename Member>
    KOKKOS_INLINE_FUNCTION static void getForces(  ExePolicy element_policy,
                                                   ExePolicy vert_policy,
                                                   T scratch,
                                                   ConstScalarElementViewType l0,
                                                   ConstScalarElementViewType l,
                                                   ConstScalarDofViewType v,
                                                   ConstScalarDofViewType current_coords,
                                                   ConnectivityViewType connectivity,
                                                   ConstScalarVertViewType mass_matrix,
                                                   ScalarDofViewType f_int,
                                                   ScalarDofViewType f,
                                                   const Member& team)
    {
      // swap force arrays and zero internal forces
      Kokkos::parallel_for(vert_policy, [=](const Ordinal i) {
        TeamOuterLoop::getForceLoop1(i, f_int);
      });
      // barrier needed between loop 1 and loop 2
      team.team_barrier();
      // load the material properties into registers so we don't
      // constantly hit the shared memory in the middle of the hot loop
      Scalar elastic_modulus = scratch(0).fiber_elastic_modulus;
      Scalar area = scratch(0).fiber_area;
      Kokkos::parallel_for(
          element_policy,
          [=](const Ordinal i) {
            TeamOuterLoop::getForceLoop2(
                i, l0, l,
                current_coords, connectivity,  f_int,elastic_modulus,area);
          });
      // barrier needed between loop 2 and loop3
      team.team_barrier();
      // load the material properties into registers so we don't
      // constantly hit the shared memory in the middle of the hot loop
      Scalar viscous_damping_coefficient = scratch(0).viscous_damping_coefficient;
      if(scratch(0).step % UPDATE_FREQ == 0)
      {
        Scalar residual = 0;
        Kokkos::parallel_reduce(
            vert_policy,
            [=](const Ordinal i, Scalar & residual_update) {
              residual_update +=
                  TeamOuterLoop::getForceLoop3(i, mass_matrix, v, f_int, f,viscous_damping_coefficient);
            },
            residual);
        Kokkos::single(Kokkos::PerTeam(team), [=]()
        {
          scratch(0).residual = sqrt(residual);
        });
      }
      else
      {
        Kokkos::parallel_for(
            vert_policy,
            [=](const Ordinal i) {
                  TeamOuterLoop::getForceLoop3(i, mass_matrix, v, f_int, f, viscous_damping_coefficient);
            });
      }
      // we don't put barrier here, because we expect the user to place abarrier after the call if needed
    }
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void updateAccelVels(ExePolicy vert_policy,
                                                       Scalar dt,
                                                       ConstScalarVertViewType mass_matrix,
                                                       ConstScalarDofViewType f,
                                                       ScalarDofViewType v)
    {
      // 2*ndof loads, 2*ndof writes
      // 2*ndof multiply, 2*ndof divide
      Kokkos::parallel_for(vert_policy, [=](const Ordinal i) {
        TeamOuterLoop::updateAccelVel(i, dt, mass_matrix, f, v);
      });
    }
    template <typename ExePolicy>
    KOKKOS_INLINE_FUNCTION static void updates(ExePolicy vert_policy,
                                               ScalarDofViewType velocity,
                                               ConstScalarVertViewType mass_matrix,
                                               ConstScalarDofViewType f,
                                               Scalar dt,
                                               ScalarDofViewType current_coordinates)
    {
      Kokkos::parallel_for(vert_policy, [=](const Ordinal i) {
         
         Scalar mass = mass_matrix(i);
         for(int j=0; j<3; ++j)
         {
           Scalar v_local = velocity(i,j) + dt * f(i,j)/mass;
           velocity(i,j) = v_local;
           current_coordinates(i,j) = current_coordinates(i,j) + 2*dt*v_local;
         }
      });
    }
    // note ever function we call has been hand tuned to reduce regsiter count so we
    // lie at the 64 reg boundary (so we can acheive 50% occupancy). Please verify
    // that you are not incraseing the reg count using the --resource-usage flag
    // when modifying this code!
    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    static bool run(Ordinal num_rves,
                    ConnectivityType connectivity,
                    T1 original_coordinates,
                    T1 current_coordinates,
                    T1 displacement,
                    T1 velocity,
                    T1 force_total,
                    T1 force_internal,
                    T2 nodal_mass,
                    T3 original_length,
                    T3 current_length,
                    T4 fiber_elastic_modulus,
                    T4 fiber_area,
                    T4 fiber_density,
                    T4 viscous_damping_coefficient,
                    T4 critical_time_scale_factor,
                    T5 displacement_boundary_vert,
                    Ordinal num_threads)
    {
      using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace, Kokkos::IndexType<Ordinal>, Kokkos::Schedule<Kokkos::Static>>;
      using member_type = typename OuterPolicyType::member_type;
      OuterPolicyType outer_policy(num_rves, num_threads);
      Kokkos::View<Scalar*,ExeSpace> dt_crit_d("dt_crit", num_rves);
      Kokkos::View<Scalar*,ExeSpace> elastic_modulus_d("elastic_modulus", num_rves);
      using ScratchPad =
          Kokkos::View<ScratchRegister<Scalar,Ordinal>*,
                       typename ExeSpace::scratch_memory_space,
                       Kokkos::MemoryUnmanaged>;
        constexpr int num_scratch = 1;
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
            auto displacement_boundary_vert_row =
                displacement_boundary_vert.template getRow<ExeSpace>(i);
            //Kokkos::single(Kokkos::PerTeam(team_member), [=]()
            // Using kokkos single uses ~3 registers here (sm_70)
            if(team_member.team_rank() == 0)
            {
                  scratch(0).viscous_damping_coefficient = 
                                viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_elastic_modulus = 
                                  fiber_elastic_modulus.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_area = 
                                  fiber_area.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_density = 
                                  fiber_density.template getRow<ExeSpace>(i)(0);
                  scratch(0).dt_crit = 0;
                  scratch(0).residual = 10.0;
                  scratch(0).step = 0;
              }
            // create the loop policies based on the lengths of the various
            // arrays
            auto element_policy = Kokkos::TeamThreadRange(
                team_member, 0, current_length_row.extent(0));
            auto vert_policy = Kokkos::TeamThreadRange(
                team_member, 0, displacement_row.extent(0));
            auto free_vert_policy = Kokkos::TeamThreadRange(
                team_member, 0,
                displacement_row.extent(0) -
                    displacement_boundary_vert_row.extent(0));
            // need to barrier here since the remaining policies need to use the scratch
            // data
            team_member.team_barrier();
            TeamOuterLoop::getElementLengths(
                element_policy, original_coordinates_row, connectivity_row,
                original_length_row, scratch, team_member);
            // very important that we update the coordinates with the initial displacements
            // the main loop only updates coordinates by using c(n+1) = c(n) + du
            TeamOuterLoop::getCurrentCoords(
                vert_policy, original_coordinates_row, displacement_row,
                current_coordinates_row);
            team_member.team_barrier();
            TeamOuterLoop::getElementLengths(
                element_policy, current_coordinates_row, connectivity_row,
                current_length_row, scratch, team_member);
            team_member.team_barrier();
            // presumably, we only need to do this over the free DOFs
            TeamOuterLoop::getForces(
                element_policy, free_vert_policy, scratch,
                original_length_row, current_length_row,
                velocity_row, current_coordinates_row, connectivity_row,
                nodal_mass_row, force_internal_row,
                force_total_row, team_member);
            team_member.team_barrier();
             dt_crit_d(i) = scratch(0).dt_crit;
          });
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
            auto displacement_boundary_vert_row =
                displacement_boundary_vert.template getRow<ExeSpace>(i);
            // we never actually change this, so save some regs and make it constexpr
            constexpr Scalar crit_time_scale_factor = 0.8;
            if(team_member.team_rank() == 0)
            {
                  scratch(0).viscous_damping_coefficient = 
                                viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_elastic_modulus = 
                                  fiber_elastic_modulus.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_area = 
                                  fiber_area.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_density = 
                                  fiber_density.template getRow<ExeSpace>(i)(0);
                  scratch(0).dt_crit = dt_crit_d(i);
                  scratch(0).residual = 10.0;
                  scratch(0).step = 0;
              }
            // create the loop policies based on the lengths of the various
            // arrays
            auto element_policy = Kokkos::TeamThreadRange(
                team_member, 0, current_length_row.extent(0));
            auto vert_policy = Kokkos::TeamThreadRange(
                team_member, 0, force_total_row.extent(0));
            auto free_vert_policy = Kokkos::TeamThreadRange(
                team_member, 0,
                force_total_row.extent(0) -
                    displacement_boundary_vert_row.extent(0));
            //  need ot have a barrier
            //  to make sure any scratch data has come in properly
            team_member.team_barrier();
            do
            {
              // Kokkos::single was causing an extra 5 registers here! (sm_70)
              if(team_member.team_rank() == 0)
              {
                scratch(0).dt_nphalf = crit_time_scale_factor * scratch(0).dt_crit*0.5;
              }
              team_member.team_barrier();
              TeamOuterLoop::updates(free_vert_policy, velocity_row, nodal_mass_row, force_total_row,
                                     scratch(0).dt_nphalf, current_coordinates_row);
              team_member.team_barrier();
              // this policy also computes the critical dt (min element length/speed of sound
              TeamOuterLoop::getElementLengths(
                  element_policy, current_coordinates_row, connectivity_row,
                  current_length_row, scratch, team_member);
              team_member.team_barrier();
              TeamOuterLoop::getForces(
                   element_policy, free_vert_policy, scratch,
                   original_length_row, current_length_row,
                  velocity_row, current_coordinates_row, connectivity_row,
                  nodal_mass_row, force_internal_row,
                  force_total_row, team_member);
              team_member.team_barrier();
              TeamOuterLoop::updateAccelVels(
                  free_vert_policy, scratch(0).dt_nphalf, nodal_mass_row,
                  force_total_row, velocity_row);
              team_member.team_barrier();
              if(team_member.team_rank() == 0)
              {
                ++scratch(0).step;
              }
            } while (scratch(0).residual > 1E-6);
          });
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
            if(team_member.team_rank() == 0)
            {
                  scratch(0).viscous_damping_coefficient = 
                                viscous_damping_coefficient.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_elastic_modulus = 
                                  fiber_elastic_modulus.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_area = 
                                  fiber_area.template getRow<ExeSpace>(i)(0);
                  scratch(0).fiber_density = 
                                  fiber_density.template getRow<ExeSpace>(i)(0);
                  scratch(0).dt_crit = dt_crit_d(i);
                  scratch(0).residual = 10.0;
                  scratch(0).step = 0;
             }
            // create the loop policies based on the lengths of the various
            // arrays
            auto element_policy = Kokkos::TeamThreadRange(
                team_member, 0, current_length_row.extent(0));
            auto vert_policy = Kokkos::TeamThreadRange(
                team_member, 0, displacement_row.extent(0));
            // barrier needed to make sure we wait for shared data
            team_member.team_barrier();
            // This extra call is required to getForces with the full vert_policy
            // to fill in the internal forces on the boundary (which is needed
            // for computation of the cauchy stresses).
              TeamOuterLoop::getForces(
                   element_policy, vert_policy, scratch,
                   original_length_row, current_length_row,
                  velocity_row, current_coordinates_row, connectivity_row,
                  nodal_mass_row, force_internal_row,
                  force_total_row, team_member);
            // Need to fill the displacement row with current_coords-initial_coords
              Kokkos::parallel_for(vert_policy, [=](const int i){
                  for(int j=0; j<3; ++j)
                  {
                    displacement_row(i,j) = current_coordinates_row(i,j) - original_coordinates_row(i,j);
                  }
              });
            team_member.team_barrier();
          });

      return true;
    }
  };
}  // namespace bio
#endif
