#include <Kokkos_Core.hpp>

namespace impl
{
  enum class BatchLoopType
  {
    SINGLE,
    TEAM_OUTER_WHILE,
    TEAM_INNER_WHILE
  };
  template<typename Scalar>
  KOKKOS_INLINE_FUNCTION
  Scalar 
  getLinearReactionForce(Scalar orig_length,
                                Scalar length,
                                Scalar elastic_modulus,
                                Scalar area)
  {
    // abaqus ...
    // double length_ratio = length / orig_length;
    // double log_strain = log(length_ratio);
    // return elastic_modulus * area * log_strain / length_ratio;
    double length_ratio = length / orig_length;
    double green_strain = 1.0 / 2.0 * (length_ratio * length_ratio - 1);
    return length_ratio * elastic_modulus * area * green_strain;
  }

  template<typename Scalar, typename LocalOrdinal, typename T1, typename T2>
  KOKKOS_INLINE_FUNCTION
  Scalar getCurrentCoordinate(LocalOrdinal node_idx,
                              LocalOrdinal direction,
                              T1 coordinates,
                              T2 displacements)
  {
    return coordinates(node_idx,direction)+displacemens(node_idx, direction);
  }

  template<typename Scalar, typename LocalOrdinal, typename T1, typename T2>
  KOKKOS_INLINE_FUNCTION
  Scalar getElementLength(LocalOrdinal element_idx,
                    T1 connectivity, 
                    T2 coordinates)
  {
    auto n1 = connectivity(element_idx, 0);
    auto n2 = connectivity(element_idx, 1);
    auto x1 = coordinates(n2,0) - coordinates(n1,0);
    auto x2 = coordinates(n2,1) - coordinates(n1,1);
    auto x3 = coordinates(n2,2) - coordinates(n1,2);
    return sqrt(x1*x1 + x2*x2 + x3*x3);
  }

  template <typename LocalOrdinal, typename T1, typename T2, typename O1>
  KOKKOS_INLINE_FUNCTION
  void updateAcceleration(LocalOrdinal node_idx, T1 mass, T2 force, O1 acceleration)
  {
    auto one_ovr_mass = 1.0/mass(node_idx);
    acceleration(node_idx,0) = one_ovr_mass * force(node_idx, 0);
    acceleration(node_idx,1) = one_ovr_mass * force(node_idx, 1);
    acceleration(node_idx,2) = one_ovr_mass * force(node_idx, 2);
  }

  template<typename LocalOrdinal, typename Scalar, typename T1, typename O1>
  KOKKOS_INLINE_FUNCTION
  void updateVelocity(LocalOrdinal node_idx, Scalar dt, T1 acceleration, O1 velocity)
  {
    velocity(node_idx,0) += dt*acceleration(node_idx, 0);
    velocity(node_idx,1) += dt*acceleration(node_idx, 1);
    velocity(node_idx,2) += dt*acceleration(node_idx, 2);
  };

  // we may want to add a displacement delta if we end up computing energies
  template<typename LocalOrdinal, typename Scalar, typename T1, typename O1>
  KOKKOS_INLINE_FUNCTION
  void updateDisplacement(LocalOrdinal node_idx, Scalar dt, T1 velocity, O1 displacement)
  {
    displacement(node_idx,0) += dt*velocity(node_idx, 0);
    displacement(node_idx,1) += dt*velocity(node_idx, 1);
    displacement(node_idx,2) += dt*velocity(node_idx, 2);
  };

  template<typename LocalOrdinal, typename Scalar, typename T1, typename O1, typename O2>
  KOKKOS_INLINE_FUNCTION
  void updateAccelerationVelocity(LocalOrdinal node_idx, Scalar dt, T1 force, O1 acceleration, O2 velocity)
  {
    auto one_ovr_mass = 1.0/mass(node_idx);
    auto a_local_1 = one_ovr_mass * force(node_idx, 0);
    auto a_local_2 = one_ovr_mass * force(node_idx, 1);
    auto a_local_3 = one_ovr_mass * force(node_idx, 2);
    acceleration(node_idx,0) = a_local_1;
    acceleration(node_idx,1) = a_local_2;
    acceleration(node_idx,2) = a_local_3;
    velocity(node_idx,0) += dt*a_local_1;
    velocity(node_idx,1) += dt*a_local_2;
    velocity(node_idx,2) += dt*a_local_3;
  };

  // Note: applying all of the boundary conditions uses the same
  // functional form, so we do not repeat these functions
  // here we do the math to back out the the node idx and direction
  // because this needs to move a smaller array to the GPU, and we are 
  // typically bandwith limited, not operation limited
  template <typename LocalOrdinal, typename Scalar, typename T1, typename T2, typename O1>
  KOKKOS_INLINE_FUNCTION
  void applyBC(LocalOrdinal fixed_idx, Scalar amplitude, T1 fixed_dofs, T2 Fixed_value, O1 update_value)
  {
    LocalOrdinal dof = fixed_dofs(fixed_idx);
    // divide by 3 assumes we are using a three dimensional mesh
    LocalOrdinal node_idx = dof/3;
    LocalOrdinal direction_idx = dof%3;

    update_value(node_idx,direction_idx) = fixed_value(fixed_idx)*amplitude;
  }

  // simultaneously updates 2 BCs...This should save one memory load
  template <typename LocalOrdinal, typename Scalar, typename T1, typename T2, typename T3, typename O1, typename O2>
  KOKKOS_INLINE_FUNCTION
  void apply2BC(LocalOrdinal fixed_idx, Scalar amplitude1, Scalar amplitude2,
                T1 fixed_dofs1, T2 fixed_value1, T3 fixed_value2,
                O1 update_value1, O2 update_value2)
  {
    LocalOrdinal dof = fixed_dofs(fixed_idx);
    // divide by 3 assumes we are using a three dimensional mesh
    LocalOrdinal node_idx = dof/3;
    LocalOrdinal direction_idx = dof%3;

    update_value1(node_idx,direction_idx) = fixed_value1(fixed_idx)*amplitude1;
    update_value2(node_idx,direction_idx) = fixed_value2(fixed_idx)*amplitude2;
  }
  // this is only needed if we care what forces are on the boundary since we
  // always have displacement BCs
  template<typename LocalOrdinal, typename T1, typename T2, typename T3, typename T4, typename T5,
    typename O1, typename O2>
  KOKKOS_INLINE_FUNCTION
  void applyBoundaryForces(LocalOrdinal fixed_idx, T1 fixed_dofs, T2 mass, T3 acceleration, T4 f_int, T5 f_damp, O1 f_ext, O2 force)
  {
    LocalOrdinal dof = fixed_dofs(fixed_idx);
    // divide by 3 assumes we are using a three dimensional mesh
    LocalOrdinal node_idx = dof/3;
    LocalOrdinal direction_idx = dof%3;

    auto f_inertial = mass(node_idx)*acceleration(node_idx,direction_idx);
    f_ext(node_idx, direction_idx) = f_inertial + f_int(node_idx,direction_idx) + f_damp(node_idx,direction_idx);
    force(node_idx, direction_idx) = f_inertial;
  }


  // TODO IMPLEMENT GET FORCES STUFF
  // call Kokkos::deep_copy(f_int, 0) before calling the main gather scatter stuff...
  
  template<typename LocalOrdinal, typename Scalar, typename T1, typename T2, typename O1>
  KOKKOS_INLINE_FUNCTION
  void updateDampingForce(LocalOrdinal node_idx, Scalar visc_damp_coeff, T1 mass, T2 velocity, O1 f_damp)
  {
    f_damp(node_idx,0) = visc_damp_coeff*mass(node_idx)*velocity(node_idx,0);
    f_damp(node_idx,1) = visc_damp_coeff*mass(node_idx)*velocity(node_idx,1);
    f_damp(node_idx,2) = visc_damp_coeff*mass(node_idx)*velocity(node_idx,2);
  }
  
  template<typename LocalOrdinal, typename Scalar, typename T1>
  KOKKOS_INLINE_FUNCTION
  Scalar getElementNorm(LocalOrdinal node1_idx, LocalOrdinal node2_idx, LocalOrdinal direction, Scalar length, T1 current_coordinates)
  {
    return (current_coordinates(node2_idx,direction) - current_coordinates(node1_idx,direction))/length;
  }

}
