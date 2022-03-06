#ifndef MUMFIM_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_SERIAL_OUTER_LOOP_H
#define MUMFIM_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_SERIAL_OUTER_LOOP_H
#include "bioBatchedFiberRVEAnalysisExplicitBase.h"

#define UPDATE_FREQ 100

namespace mumfim
{
  template <typename Scalar, typename Ordinal, typename ExeSpace>
  struct SerialOuterLoop
      : public BaseBatchedExplicit<SerialOuterLoop<Scalar, Ordinal, ExeSpace>,
                                   Scalar,
                                   Ordinal,
                                   ExeSpace>
  {
    using BaseType =
        BaseBatchedExplicit<SerialOuterLoop<Scalar, Ordinal, ExeSpace>,
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
    static void getCurrentCoords(ExePolicy vert_policy,
                                 ConstScalarDofViewType coords,
                                 ConstScalarDofViewType u,
                                 ScalarDofViewType current_coords)
    {
      Kokkos::parallel_for(
          "getCurrentCoords", vert_policy, KOKKOS_LAMBDA(const Ordinal i) {
            SerialOuterLoop::getCurrentCoord(i, coords, u, current_coords);
          });
    }
    template <typename ExePolicy>
    static void getElementLengths(ExePolicy element_policy,
                                  ConstScalarDofViewType coords,
                                  ConnectivityViewType connectivity,
                                  ScalarElementViewType l0,
                                  Ordinal step,
                                  Scalar fiber_elastic_modulus,
                                  Scalar fiber_density,
                                  Scalar& dt_crit)
    {
      if(step % UPDATE_FREQ == 0)
      {
        Scalar l_min;
        Kokkos::Min<Scalar> min_reducer(l_min);
        Kokkos::parallel_reduce(element_policy, KOKKOS_LAMBDA(Ordinal i, Scalar& local) {
          SerialOuterLoop::getElementLength(i, coords, connectivity, l0);
          min_reducer.join(local, l0(i));
        }, min_reducer);
        Scalar sound_speed = sqrt(fiber_elastic_modulus / fiber_density);
        dt_crit = l_min/sound_speed;
      }
      else
      {
        Kokkos::parallel_for(element_policy, KOKKOS_LAMBDA(Ordinal i) {
          SerialOuterLoop::getElementLength(i, coords, connectivity, l0);
        });
      }
    }
    template <typename ExePolicy>
    static void updateAccelVels(ExePolicy vert_policy,
                                Scalar dt,
                                ConstScalarVertViewType mass_matrix,
                                ConstScalarDofViewType f,
                                ScalarDofViewType v)
    {
      // 2*ndof loads, 2*ndof writes
      // 2*ndof multiply, 2*ndof divide
      Kokkos::parallel_for(
          "updateAccelVel", vert_policy, KOKKOS_LAMBDA(const int i) {
            SerialOuterLoop::updateAccelVel(i, dt, mass_matrix, f, v);
          });
    }
    template <typename ExePolicy>
    static void updates(ExePolicy vert_policy,
                                               ScalarDofViewType velocity,
                                               ConstScalarVertViewType mass_matrix,
                                               ConstScalarDofViewType f,
                                               Scalar dt,
                                               ScalarDofViewType current_coordinates)
    {
      Kokkos::parallel_for(vert_policy, KOKKOS_LAMBDA(const Ordinal i) {
         
         Scalar mass = mass_matrix(i);
         for(int j=0; j<3; ++j)
         {
           Scalar v_local = velocity(i,j) + dt * f(i,j)/mass;
           velocity(i,j) = v_local;
           current_coordinates(i,j) = current_coordinates(i,j) + 2*dt*v_local;
         }
      });
    }
    template <typename ExePolicy>
    static void getForces(  ExePolicy element_policy,
                                                   ExePolicy vert_policy,
                                                   ConstScalarElementViewType l0,
                                                   ConstScalarElementViewType l,
                                                   ConstScalarDofViewType v,
                                                   ConstScalarDofViewType current_coords,
                                                   ConnectivityViewType connectivity,
                                                   ConstScalarVertViewType mass_matrix,
                                                   ScalarDofViewType f_int,
                                                   ScalarDofViewType f,
                                                   Scalar fiber_elastic_modulus,
                                                   Scalar fiber_area,
                                                   Scalar viscous_damping_coefficient,
                                                   Ordinal step,
                                                   Scalar & residual)
    {
      // swap force arrays and zero internal forces
      Kokkos::parallel_for(vert_policy, KOKKOS_LAMBDA(const Ordinal i) {
        SerialOuterLoop::getForceLoop1(i, f_int);
      });
      Kokkos::parallel_for(
          element_policy,
          KOKKOS_LAMBDA(const Ordinal i) {
            SerialOuterLoop::getForceLoop2(
                i, l0, l,
                current_coords, connectivity,  f_int,fiber_elastic_modulus,fiber_area);
          });
      if(step % UPDATE_FREQ == 0)
      {
        residual = 0;
        Kokkos::parallel_reduce(
            vert_policy,
            KOKKOS_LAMBDA(const Ordinal i, Scalar & residual_update) {
              residual_update +=
                  SerialOuterLoop::getForceLoop3(i, mass_matrix, v, f_int, f,viscous_damping_coefficient);
            },
            residual);
          residual = sqrt(residual);
      }
      else
      {
        Kokkos::parallel_for(
            vert_policy,
            KOKKOS_LAMBDA(const Ordinal i) {
                  SerialOuterLoop::getForceLoop3(i, mass_matrix, v, f_int, f, viscous_damping_coefficient);
            });
      }
    }
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
      for (size_t i = 0; i < num_rves; ++i)
      {
        //// get the device views of the data
        auto connectivity_row = connectivity.template getRow<ExeSpace>(i);
        auto original_coordinates_row =
            original_coordinates.template getRow<ExeSpace>(i);
        auto current_coordinates_row =
            current_coordinates.template getRow<ExeSpace>(i);
        auto displacement_row = displacement.template getRow<ExeSpace>(i);
        auto velocity_row = velocity.template getRow<ExeSpace>(i);
        auto force_total_row = force_total.template getRow<ExeSpace>(i);
        auto force_internal_row = force_internal.template getRow<ExeSpace>(i);
        auto nodal_mass_row = nodal_mass.template getRow<ExeSpace>(i);
        auto original_length_row = original_length.template getRow<ExeSpace>(i);
        auto current_length_row = current_length.template getRow<ExeSpace>(i);
        auto displacement_boundary_vert_row =
            displacement_boundary_vert.template getRow<ExeSpace>(i);
        constexpr Scalar crit_time_scale_factor = 0.8;
        Scalar elastic_modulus =
            fiber_elastic_modulus.template getRow<Kokkos::HostSpace>(i)(0);
        Scalar area = fiber_area.template getRow<Kokkos::HostSpace>(i)(0);
        Scalar density = fiber_density.template getRow<Kokkos::HostSpace>(i)(0);
        Scalar visc_damp_coeff = viscous_damping_coefficient.template getRow<Kokkos::HostSpace>(i)(0);
        // create the loop policies based on the lengths of the various arrays
        Kokkos::RangePolicy<ExeSpace> element_policy(
            0, current_length_row.extent(0));
        Kokkos::RangePolicy<ExeSpace> vert_policy(0, displacement_row.extent(0));
        Kokkos::RangePolicy<ExeSpace> free_vert_policy(0,
            displacement_row.extent(0) -
                displacement_boundary_vert_row.extent(0));
        Scalar residual = 10.0;
        Scalar dt_nphalf;
        long int n_step = 0;
        Scalar dt_crit;
        SerialOuterLoop::getElementLengths(
            element_policy, original_coordinates_row, connectivity_row,
            original_length_row, n_step, elastic_modulus, density, dt_crit);
        SerialOuterLoop::getCurrentCoords(vert_policy, original_coordinates_row,
                                          displacement_row,
                                          current_coordinates_row);
        SerialOuterLoop::getElementLengths(
            element_policy, current_coordinates_row, connectivity_row,
            current_length_row,n_step, elastic_modulus, density, dt_crit);

          SerialOuterLoop::getForces(
               element_policy, free_vert_policy,
               original_length_row, current_length_row,
              velocity_row, current_coordinates_row, connectivity_row,
              nodal_mass_row, force_internal_row,
              force_total_row,
              elastic_modulus, area, visc_damp_coeff, n_step,residual);
        do
        {
          dt_nphalf = crit_time_scale_factor * dt_crit*0.5;
          SerialOuterLoop::updates(free_vert_policy, velocity_row, nodal_mass_row,
                                   force_total_row,
                                   dt_nphalf, current_coordinates_row);
          SerialOuterLoop::getElementLengths(
              element_policy, current_coordinates_row, connectivity_row,
              current_length_row, n_step, elastic_modulus, density, dt_crit);
          SerialOuterLoop::getForces(
               element_policy, free_vert_policy,
               original_length_row, current_length_row,
              velocity_row, current_coordinates_row, connectivity_row,
              nodal_mass_row, force_internal_row,
              force_total_row,
              elastic_modulus, area, visc_damp_coeff, n_step,residual);
          SerialOuterLoop::updateAccelVels(
              free_vert_policy, dt_nphalf, nodal_mass_row,
              force_total_row, velocity_row);
          ++n_step;
          Kokkos::fence();
        } while (residual > 1E-6);
      }
      return true;
    }
  };
}  // namespace mumfim
#endif
