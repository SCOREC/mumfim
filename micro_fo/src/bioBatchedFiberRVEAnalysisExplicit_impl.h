#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_IMPL_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_IMPL_H__
#include <apf.h>           // for extractCoordinateArray
#include <apfConvert.h>    // for extractCoordinateArray
#include <apfMesh2.h>      // for extractCoordinateArray
#include <apfNumbering.h>  // for getNumbering (getFixedDOF)
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <array>
#include <iostream>
#include <limits>
#include "bioMassIntegrator.h"
#include "bioMicroFOParams.h"
#include "bioRVE.h"
#include "bioUtility.h"
namespace bio
{
  template <typename T>
  void getConnectivity(apf::Mesh * m, T connectivity_view, int cellDim)
  {
    bool first = true;
    auto nelem = m->count(cellDim);
    apf::GlobalNumbering * global =
        apf::makeGlobal(apf::numberOwnedNodes(m, "apf_destruct"));
    apf::synchronize(global);
    apf::MeshIterator * it = m->begin(cellDim);
    apf::MeshEntity * e;
    int i = 0;
    while ((e = m->iterate(it)))
    {
      auto etype = m->getType(e);
      apf::Downward verts;
      int nverts = m->getDownward(e, 0, verts);
      if (first)
      {
        first = false;
        if (connectivity_view.extent(0) != nelem * nverts)
        {
          std::cerr << "The connectivity_view has "
                    << connectivity_view.extent(0)
                    << "elements, but it should have " << nelem * nverts
                    << "elements." << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      for (int j = 0; j < nverts; ++j)
      {
        connectivity_view(i++) = apf::getNumber(global, apf::Node(verts[j], 0));
      }
    }
    m->end(it);
    apf::destroyGlobalNumbering(global);
  }
  template <typename T, int NumValues = 3>
  void getFieldValues(apf::Mesh * m, apf::Field * field, T field_view)
  {
    apf::MeshIterator * it = m->begin(0);
    std::array<double, NumValues> values;
    int i = 0;
    while (apf::MeshEntity * v = m->iterate(it))
    {
      if (m->isOwned(v))
      {
        apf::getComponents(field, v, 0, values.data());
        for (int v = 0; v < NumValues; ++v)
        {
          field_view(i * NumValues + v) = values[v];
        }
        i++;
      }
    }
    m->end(it);
  }
  // TODO this can be made more efficient by having the integrator load
  // direcetly into the Kokkos::View, but we will hold off on this
  // because ultimately we need device side integration routines
  template <typename T>
  void getNodalMass(apf::Mesh * mesh,
                    double fiber_density,
                    double fiber_area,
                    T mass_matrix_view)
  {
    auto nodal_mass =
        apf::createLagrangeField(mesh, "nodalMass_tmp", apf::SCALAR, 1);
    apf::zeroField(nodal_mass);
    // build the mass integrator and load the nodal mass into a temporary field
    MassIntegrator massInt(nodal_mass, fiber_density, fiber_area, 3,
                           MassLumpType::RowSum);
    massInt.process(mesh, 1);
    // copy the field values to the kokkos view
    getFieldValues<T, 1>(mesh, nodal_mass, mass_matrix_view);
    apf::destroyField(nodal_mass);
  }
  template <typename T>
  void getFixedDof(apf::Numbering * numbering,
                   T fixed_dof_view,
                   std::vector<apf::MeshEntity *> bnd_nds)
  {
    for (std::size_t i = 0; i < bnd_nds.size(); ++i)
    {
      apf::MeshEntity * nd = bnd_nds[i];
      fixed_dof_view(i * 3) = apf::getNumber(numbering, nd, 0, 0);
      fixed_dof_view(i * 3 + 1) = apf::getNumber(numbering, nd, 0, 1);
      fixed_dof_view(i * 3 + 2) = apf::getNumber(numbering, nd, 0, 2);
    }
  }
  template <typename T>
  void updateRVECoords(RVE & rve, T incremental_deformation_gradient)
  {
    apf::Matrix3x3 F;
    for (int ei = 0; ei < 3; ++ei)
      for (int ej = 0; ej < 3; ++ej)
        F[ei][ej] = incremental_deformation_gradient(ei, ej);
    ApplyIncrementalDeformationGradient(F, rve.getMesh(), rve.getdUField(),
                                        rve.getUField());
  }
  template<typename ExeSpace, typename T, typename T2, typename T3, int DIM=3>
  void applyIncrementalDeformationToDisplacement(T deformation_gradient, T3 coordinates, T2 displacements)
  {
    deformation_gradient.template sync<ExeSpace>();
    auto deformation_gradient_d = deformation_gradient.template view<ExeSpace>();
    coordinates.template sync<ExeSpace>();
    displacements.template sync<ExeSpace>();
    displacements.template modify<ExeSpace>();
    using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace>;
    using member_type = typename OuterPolicyType::member_type;
    OuterPolicyType outer_policy(deformation_gradient.extent(0), Kokkos::AUTO());

    Kokkos::parallel_for("applyIncrementalDfm", outer_policy, KOKKOS_LAMBDA(member_type team_member){
        int rve = team_member.league_rank();
        auto displacements_row = displacements.template getRow<ExeSpace>(rve);
        auto coordinates_row = coordinates.template getRow<ExeSpace>(rve);
        auto deformation_gradient_row = Kokkos::subview(deformation_gradient_d, rve, Kokkos::ALL(), Kokkos::ALL());
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, displacements_row.extent(0)/DIM),[&](int i) {
            for(int k=0; k<3; ++k)
            {
            displacements_row(i*DIM+k) = (deformation_gradient_row(k,0)-(k==0?1:0))*coordinates_row(i*DIM) +
                                       (deformation_gradient_row(k,1)-(k==1?1:0))*coordinates_row(i*DIM+1) +
                                       (deformation_gradient_row(k,2)-(k==2?1:0))*coordinates_row(i*DIM+2) +
                                       displacements_row(i*DIM+k);
            }
           });
        });

  }
  template<typename ExeSpace, typename T, typename T2, typename T3, int DIM=3>
  void applyIncrementalDeformationToBoundary(T deformation_gradient, T2 coordinates,
                                             T3 boundary_dof, T2 boundary_values)
  {
    deformation_gradient.template sync<ExeSpace>();
    auto deformation_gradient_d = deformation_gradient.template view<ExeSpace>();
    coordinates.template sync<ExeSpace>();
    boundary_dof.template sync<ExeSpace>();
    boundary_values.template sync<ExeSpace>();
    boundary_values.template modify<ExeSpace>();

    using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace>;
    using member_type = typename OuterPolicyType::member_type;
    OuterPolicyType outer_policy(deformation_gradient.extent(0), Kokkos::AUTO());

    Kokkos::parallel_for("applyIncrementalDfm", outer_policy, KOKKOS_LAMBDA(member_type team_member){
        int rve = team_member.league_rank();
        auto deformation_gradient_row = Kokkos::subview(deformation_gradient_d, rve, Kokkos::ALL(), Kokkos::ALL());
        auto coordinates_row = coordinates.template getRow<ExeSpace>(rve);
        auto boundary_dof_row = boundary_dof.template getRow<ExeSpace>(rve);
        auto boundary_values_row = boundary_values.template getRow<ExeSpace>(rve);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, boundary_dof_row.extent(0)/DIM),[&](int i) {
            auto dof1 = boundary_dof_row(i*DIM);
            auto dof2 = boundary_dof_row(i*DIM+1);
            auto dof3 = boundary_dof_row(i*DIM+2);
            for(int k=0; k<3; ++k)
            {
            boundary_values_row(i*DIM+k) = (deformation_gradient_row(k,0)-(k==0?1:0))*coordinates_row(dof1) +
                                       (deformation_gradient_row(k,1)-(k==1?1:0))*coordinates_row(dof2) +
                                       (deformation_gradient_row(k,2)-(k==2?1:0))*coordinates_row(dof3) +
                                       boundary_values_row(i*DIM+k);
            }
           });
        });

  }
  template<typename ExeSpace, typename T, typename T2, typename T3>
  void computeCauchyStress(T2 boundary_dofs, T3 coordinates, T3 force, T stress)
  {
    boundary_dofs.template sync<ExeSpace>();
    coordinates.template sync<ExeSpace>();
    force.template sync<ExeSpace>();

    stress.template modify<ExeSpace>();
    auto stress_d = stress.template view<ExeSpace>();
    // zero the trss vector because we will sum the values into the streses
    Kokkos::deep_copy(stress_d, 0);

    using OuterPolicyType = Kokkos::TeamPolicy<ExeSpace>;
    using member_type = typename OuterPolicyType::member_type;
    OuterPolicyType outer_policy(stress.extent(0), Kokkos::AUTO());

    Kokkos::parallel_for("compute cauchy stress", outer_policy, KOKKOS_LAMBDA(member_type team_member){
        int rve = team_member.league_rank();
        auto stress_row = Kokkos::subview(stress_d, rve, Kokkos::ALL());
        auto boundary_dof_row = boundary_dofs.template getRow<ExeSpace>(rve);
        auto coordinates_row = coordinates.template getRow<ExeSpace>(rve);
        auto force_row = force.template getRow<ExeSpace>(rve);

        auto loop_size = boundary_dof_row.extent(0)/3;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, loop_size),[=](const int i) {
            auto dof1 = boundary_dof_row(i*3);
            auto dof2 = boundary_dof_row(i*3+1);
            auto dof3 = boundary_dof_row(i*3+2);
            auto f1 = force_row(dof1);
            auto f2 = force_row(dof2);
            auto f3 = force_row(dof3);
            auto crd1 = coordinates_row(dof1);
            auto crd2 = coordinates_row(dof2);
            auto crd3 = coordinates_row(dof3);
            Kokkos::atomic_fetch_add(&stress_row(0), crd1*f1);
            Kokkos::atomic_fetch_add(&stress_row(1), crd2*f2);
            Kokkos::atomic_fetch_add(&stress_row(2), crd3*f3);
            Kokkos::atomic_fetch_add(&stress_row(3), 0.5*(crd2*f3+crd3*f2));
            Kokkos::atomic_fetch_add(&stress_row(4), 0.5*(crd1*f3+crd3*f1));
            Kokkos::atomic_fetch_add(&stress_row(5), 0.5*(crd1*f2+crd2*f1));
           });
        });

  }

    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar getLinearReactionForce(Scalar orig_length,
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
  template <typename ExePolicy, typename RWSV>
  void zeroDeviceData(ExePolicy exe_policy, RWSV data)
  {
    Kokkos::parallel_for(
        "zeroDeviceData", exe_policy,
        KOKKOS_LAMBDA(const int i) { data(i) = 0; });
  }
  template <typename ExePolicy, typename ROSV, typename RWSV>
  void getCurrentCoords(ExePolicy dof_policy,
                        ROSV coords,
                        ROSV u,
                        RWSV current_coords)
  {
    Kokkos::parallel_for(
        "getCurrentCoords", dof_policy,
        KOKKOS_LAMBDA(const int i) { current_coords(i) = coords(i) + u(i); });
  }
  template <typename ExePolicy, typename ROSV, typename ROOV, typename RWSV>
  void getElementLengths(ExePolicy element_policy,
                         ROSV coords,
                         ROOV connectivity,
                         RWSV l0)
  {
    Kokkos::parallel_for(
        "getElementLengths", element_policy, KOKKOS_LAMBDA(const int i) {
          int n1 = connectivity(i * 2);
          int n2 = connectivity(i * 2 + 1);
          double x1 = coords(n2 * 3) - coords(n1 * 3);
          double x2 = coords(n2 * 3 + 1) - coords(n1 * 3 + 1);
          double x3 = coords(n2 * 3 + 2) - coords(n1 * 3 + 2);
          l0(i) = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
        });
  }
  template <typename ExePolicy,
            typename Scalar,
            typename ROSV,
            typename ROOV,
            typename RWSV>
  double getForces(ExePolicy element_policy,
                   ExePolicy dof_policy,
                   Scalar fiber_elastic_modulus,
                   Scalar fiber_area,
                   Scalar fiber_density,
                   ROSV l0,
                   ROSV l,
                   ROSV v,
                   ROSV current_coords,
                   ROOV connectivity,
                   ROSV mass_matrix,
                   Scalar visc_damp_coeff,
                   RWSV f_int,
                   RWSV f_ext,
                   RWSV f_damp,
                   RWSV f,
                   Scalar & residual)
  {
    Scalar dt = std::numeric_limits<Scalar>::max();
    // element normal vectors
    Scalar sound_speed = sqrt(fiber_elastic_modulus / fiber_density);
    residual = 0;
    // swap force arrays and zero internal forces
    Kokkos::parallel_for(
        "getFoces--Loop1", dof_policy, KOKKOS_LAMBDA(const int i) {
          f_int(i) = 0;
          f_damp(i) = visc_damp_coeff * mass_matrix(i / 3) * v(i);
        });
    // set the internal forces
    Kokkos::Min<Scalar> min_reducer(dt);
    Kokkos::parallel_reduce(
        "getForces--mainLoop", element_policy,
        KOKKOS_LAMBDA(const int i, Scalar & dt_crit_elem) {
          Scalar local_l = l(i);
          // FIXME the force reactions are all messed up! need to figure out
          // how to get the struct data in here?
          Scalar frc = getLinearReactionForce(
              l0(i), local_l, fiber_elastic_modulus, fiber_area);
          auto n1 = connectivity(i * 2);
          auto n2 = connectivity(i * 2 + 1);
          Scalar elem_nrm_1 =
              (current_coords(n2 * 3) - current_coords(n1 * 3)) / local_l;
          Scalar elem_nrm_2 =
              (current_coords(n2 * 3 + 1) - current_coords(n1 * 3 + 1)) /
              local_l;
          Scalar elem_nrm_3 =
              (current_coords(n2 * 3 + 2) - current_coords(n1 * 3 + 2)) /
              local_l;
          // note we have a race condition here unless we perform an atomic
          // add!
          Kokkos::atomic_add(&f_int(n1 * 3), -frc * elem_nrm_1);
          Kokkos::atomic_add(&f_int(n1 * 3 + 1), -frc * elem_nrm_2);
          Kokkos::atomic_add(&f_int(n1 * 3 + 2), -frc * elem_nrm_3);
          Kokkos::atomic_add(&f_int(n2 * 3), frc * elem_nrm_1);
          Kokkos::atomic_add(&f_int(n2 * 3 + 1), frc * elem_nrm_2);
          Kokkos::atomic_add(&f_int(n2 * 3 + 2), frc * elem_nrm_3);
          min_reducer.join(dt_crit_elem, local_l / sound_speed);
        },
        min_reducer);
    residual = 0;
    Kokkos::parallel_reduce(
        "getForces-Loop3", dof_policy,
        KOKKOS_LAMBDA(const int i, Scalar & residual_update) {
          double local_residual = f_ext(i) - f_int(i);
          residual_update += local_residual * local_residual;
          f(i) = local_residual - f_damp(i);
        },
        residual);
    residual = sqrt(residual);
    return dt;
  }
  template <typename ExePolicy, typename ROSV, typename RWSV>
  void updateAcceleration(ExePolicy dof_policy,
                          ROSV mass_matrix,
                          ROSV f,
                          RWSV a)
  {
    Kokkos::parallel_for(
        "updateAcceleration", dof_policy, KOKKOS_LAMBDA(const int i) {
          a(i) = (1.0 / mass_matrix(i / 3)) * f(i);
        });
  }
  template <typename ExePolicy, typename ROSV, typename RWSV>
  void updateVelocity(ExePolicy dof_policy, ROSV a, double dt, RWSV v)
  {
    Kokkos::parallel_for(
        "updateVelocity", dof_policy,
        KOKKOS_LAMBDA(const int i) { v(i) += dt * a(i); });
  }
  template <typename ExePolicy, typename Scalar, typename ROSV, typename RWSV>
  void updateAccelVel(ExePolicy dof_policy,
                      Scalar dt,
                      ROSV mass_matrix,
                      ROSV f,
                      RWSV a,
                      RWSV v)
  {
    // 2*ndof loads, 2*ndof writes
    // 2*ndof multiply, 2*ndof divide
    Kokkos::parallel_for(
        "updateAccelVel", dof_policy, KOKKOS_LAMBDA(const int i) {
          Scalar a_local = (1.0 / mass_matrix(i / 3)) * f(i);
          a(i) = a_local;
          v(i) += dt * a_local;
        });
  }
  template <typename ExePolicy, typename Scalar, typename ROSV, typename RWSV>
  void updateDisplacement(ExePolicy dof_policy, ROSV v, Scalar dt, RWSV u)
  {
    Kokkos::parallel_for(
        "updateDisplacement", dof_policy,
        KOKKOS_LAMBDA(const int i) { u(i) += dt * v(i); });
  }
  template <typename ExePolicy,
            typename Scalar,
            typename ROSV,
            typename ROOV,
            typename RWSV>
  // FIXME Add init values back!
  void applyDispBC(ExePolicy fixed_dof_policy,
                   Scalar amp_t,
                   ROOV dof,
                   ROSV values,
                   RWSV u)
  {
    Kokkos::parallel_for(
        "applyDispBC", fixed_dof_policy, KOKKOS_LAMBDA(const int i) {
          u(dof(i)) = values(i) * amp_t;
        });
  }
  template <typename ExePolicy,
            typename Scalar,
            typename ROSV,
            typename ROOV,
            typename RWSV>
  void applyVelBC(ExePolicy fixed_dof_policy,
                  Scalar amp_t,
                  ROOV dof,
                  ROSV values,
                  RWSV v)
  {
    Kokkos::parallel_for(
        "applyVelBC", fixed_dof_policy,
        KOKKOS_LAMBDA(const int i) { v(dof(i)) = values(i) * amp_t; });
  }
  template <typename ExePolicy,
            typename Scalar,
            typename ROSV,
            typename ROOV,
            typename RWSV>
  void applyAccelBC(ExePolicy fixed_dof_policy,
                    Scalar amp_t,
                    ROOV dof,
                    ROSV values,
                    RWSV a)
  {
    Kokkos::parallel_for(
        "applyAccelBC", fixed_dof_policy,
        KOKKOS_LAMBDA(const int i) { a(dof(i)) = values(i) * amp_t; });
  }
  template <typename ExePolicy,
            typename Scalar,
            typename ROSV,
            typename ROOV,
            typename RWSV>
  void applyAccelVelBC(ExePolicy fixed_dof_policy,
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
    Kokkos::parallel_for(
        "applyAccelVelBC", fixed_dof_policy, KOKKOS_LAMBDA(const int i) {
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
        });
  }
}  // namespace bio
#endif
