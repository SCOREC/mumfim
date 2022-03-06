// This file is derived from ETFEM which can be found at
// https://github.com/jacobmerson/ETFEM/ The file is seperated out, so we can
// track changes in the ETFEM repository
#include <PCU.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "bioExplicitAmplitude.h"
#include "bioExplicitOutputWriter.h"
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysisExplicit.h"
#include "bioMassIntegrator.h"
#include "bioVerbosity.h"
#include <Kokkos_Core.hpp>
namespace bio
{
  KOKKOS_INLINE_FUNCTION
  static bool isClose(double a,
                      double b,
                      double rtol = 1E-8,
                      double atol = 1E-10)
  {
    return fabs(a - b) <= fmax(rtol * fmax(fabs(a), fabs(b)), atol);
  }
  static void extractCoordinateArray(apf::Mesh * m,
                                     apf::Field * coord_field,
                                     double *& coords,
                                     int & nverts)
  {
    nverts = apf::countOwned(m, 0);
    coords = new double[nverts * 3];
    apf::MeshIterator * it = m->begin(0);
    int i = 0;
    while (apf::MeshEntity * v = m->iterate(it))
    {
      if (m->isOwned(v))
      {
        // Vector3 p;
        // m->getPoint(v, 0, p);
        // p.toArray(&coords[i * 3]);
        apf::getComponents(coord_field, v, 0, &coords[i * 3]);
        i++;
      }
    }
    m->end(it);
  }
  // types
  // H_RODA: read-only double array
  // H_RWDA: read-write double array
  // H_ROIA: read-only int array
  // H_RWIA: read-write int array
  // anything with RA prefix means random access
  template <typename T,
            typename H_RODA,
            typename H_RWDA,
            typename H_ROIA,
            typename H_RWIA,
            typename D_RODA,
            typename D_RWDA,
            typename D_ROIA,
            typename D_RWIA,
            typename D_RARODA,
            typename D_RAROIA>
  class ExplicitAnalysisBase
  {
    public:
    bool run(unsigned long & itr_prev)
    {
      double residual = 1.0;
      // optimization to only compute the mass matrix 1
      // if(!mass_field_initialized)
      computeMassMatrix(nnds, mesh, nodalMass, fiber_density, fiber_area);
      // copy other data to device
      copyData(mass_matrix_d, mass_matrix_h);
      copyData(coords_d, coords_h);
      copyData(connectivity_d, connectivity_h);
      // copy over the initial velocities, accelerations, and displacements
      // this can probably be optimized out eventually
      copyData(u_d, u_h);
      copyData(v_d, v_h);
      copyData(a_d, a_h);
      copyData(f_int_d, f_int_h);
      copyData(f_ext_d, f_ext_h);
      // copy boundary data over to the device
      disp_boundary_dof_d = createDeviceIntMirrorArray(disp_boundary_dof_h);
      disp_boundary_values_d =
          createDeviceDoubleMirrorArray(disp_boundary_values_h);
      disp_boundary_init_values_d =
          createDeviceDoubleMirrorArray(disp_boundary_init_values_h);
      copyData(disp_boundary_dof_d, disp_boundary_dof_h);
      copyData(disp_boundary_values_d, disp_boundary_values_h);
      copyData(disp_boundary_init_values_d, disp_boundary_init_values_h);
      if (disp_nfixed <= 0)
      {
        std::cerr
            << "Currently we are only set up to compute problems with "
               "fixed displacement boundary conditions. Please apply some "
               "displacements.\n";
        return false;
      }
      BIO_V3(std::cout << "This problem has: " << nnds << " nodes and " << nelem
                       << " elements\n";
             std::cout << "Using Deformation Gradient BC has fixed "
                       << disp_nfixed << " of " << ndof
                       << " degrees of freedom.\n";)
      getElementLengths(nelem, coords_d, connectivity_d, l0_d);
      getCurrentCoords(ndof, coords_d, u_d, current_coords_d);
      getElementLengths(nelem, current_coords_d, connectivity_d, l_d);
      dt_crit = getForces(nelem, ndof, l0_d, l_d, v_d, current_coords_d,
                          connectivity_d, mass_matrix_d, visc_damp_coeff,
                          f_int_d, f_int_last_d, f_ext_d, f_ext_last_d,
                          f_damp_d, f_damp_last_d, f_d, residual);
      // else
      BIO_V3(std::cout << "The initial timestep is dt=" << dt_crit
                       << ".\nIf the analysis continues with a comprable dt "
                          "the analysis "
                          "will take approximately "
                       << (unsigned long)(total_time /
                                          (dt_crit * crit_time_scale_factor))
                       << " timesteps.\n";)
      // don't compute the frame time if the frequency is zero because this will
      // be a divide by zero
      double frame_time =
          print_field_frequency ? (total_time / print_field_frequency) : 0;
      unsigned int frame_num = 1;
      if (print_field_frequency)
      {
        writer->writeFieldData(itr_prev, current_time);
      }
      if (print_history_frequency)
      {
        writer->writeHistoryData(itr_prev, W_int, W_ext, W_damp, W_kin, current_time);
      }
      updateAcceleration(ndof, mass_matrix_d, f_d, a_d);

      // we set the current time to the total time
      // because we don't want to ramp boundaries.
      // The hope is that with an affine guess this
      // won't make the solve explode...
      current_time = total_time;
      do
      {
        dt_nphalf = crit_time_scale_factor * dt_crit;
        t_npone = current_time + dt_nphalf;
        // make sure the last time step is at the requested time
        //if (t_npone > total_time)
        //{
        //  t_npone = total_time;
        //  dt_nphalf = total_time - current_time;
        //}
        t_nphalf = 0.5 * (current_time + t_npone);
        updateVelocity(ndof, a_d, (t_nphalf - current_time), v_d);
        // we don't really need to apply the velocity boundary condition
        // here since the only place the boundary velocities are used are
        // to compute the damping force on the boundary which doesn't matter
        // since we are applying displacement boundary conditions.
        //applyVelBC(disp_nfixed, t_npone, disp_boundary_dof_d,
        //           disp_boundary_values_d, amp, v_d);
        updateDisplacement(ndof, v_d, dt_nphalf, u_d, du_d);
        applyDispBC(disp_nfixed, current_time, t_npone, disp_boundary_dof_d,
                    disp_boundary_init_values_d, disp_boundary_values_d, amp,
                    u_d, du_d);
        getCurrentCoords(ndof, coords_d, u_d, current_coords_d);
        getElementLengths(nelem, current_coords_d, connectivity_d, l_d);
        dt_crit = getForces(nelem, ndof, l0_d, l_d, v_d, current_coords_d,
                            connectivity_d, mass_matrix_d, visc_damp_coeff,
                            f_int_d, f_int_last_d, f_ext_d, f_ext_last_d,
                            f_damp_d, f_damp_last_d, f_d, residual);
        updateAccelVel(ndof, t_npone-t_nphalf, mass_matrix_d, f_d, a_d, v_d);
        applyAccelVelBC(disp_nfixed, t_npone, visc_damp_coeff, disp_boundary_dof_d,
                        disp_boundary_values_d, mass_matrix_d, amp, a_d, v_d,
                        f_int_d, f_ext_d, f_damp_d, f_d);
        computeEnergies(ndof, du_d, f_int_last_d, f_int_d, f_ext_last_d,
                        f_ext_d, f_damp_last_d, f_damp_d, W_int, W_ext, W_damp);
        W_kin = computeKineticEnergy(ndof, mass_matrix_d, v_d);
        // make sure we update the time before we write data to file
        ++n_step;
        current_time = t_npone;
        // make sure that all of the kernels have finished before copying data
        // back to host
        fence();
        // check balance of energy
        if (n_step > 10 &&
            !checkEnergyBalance(W_kin, W_int, W_ext, energy_check_eps))
        {
          std::cerr << "Energy Balance Not conserved in step " << n_step
                    << ". aborting analysis.\n";
          // need to copy any fields that the writer will use back to the host.
          // note that we don't need to update the mass matrix, coords, or
          // connectivity here since they do not change during the simulation
          copyData(u_h, u_d);
          copyData(v_h, v_d);
          copyData(a_h, a_d);
          copyData(f_int_h, f_int_d);
          copyData(f_ext_d, f_ext_h);
          writer->writeFieldData(itr_prev + n_step, current_time);
          writer->writeHistoryData(itr_prev + n_step, W_int, W_ext, W_damp, W_kin, current_time);
          return false;
        }
        bool print_field =
            print_field_frequency &&
            (print_field_by_num_frames
                 ? ((current_time - dt_nphalf) < frame_num * frame_time) &&
                       ((current_time >= frame_num * frame_time) ||
                        current_time >= total_time)
                 : n_step % print_field_frequency == 0);
        if (print_field)
        {
          // need to copy any fields that the writer will use back to the host.
          // note that we don't need to update the mass matrix, coords, or
          // connectivity here since they do not change during the simulation
          copyData(u_h, u_d);
          copyData(v_h, v_d);
          copyData(a_h, a_d);
          copyData(f_int_h, f_int_d);
          copyData(f_ext_d, f_ext_h);
          writer->writeFieldData(itr_prev + n_step, current_time);
          ++frame_num;
        }
        if (print_history_frequency && (n_step % print_history_frequency == 0))
        {
          writer->writeHistoryData(itr_prev + n_step, W_int, W_ext, W_damp, W_kin, current_time);
        }
        //if((current_time > total_time) && (n_step%1000 == 0))
        //  std::cout<<residual<<std::endl;
      } while ((current_time < total_time) || (residual > 1E-6));
      // if we didn't already write the last frame
      if (print_field_frequency && !print_field_by_num_frames &&
          (n_step % print_field_frequency))
      {
        copyData(u_h, u_d);
        copyData(v_h, v_d);
        copyData(a_h, a_d);
        copyData(f_int_h, f_int_d);
        copyData(f_ext_d, f_ext_h);
        writer->writeFieldData(itr_prev + n_step, current_time);
      }
      // make sure we copy the data out to the mesh even
      // if we don't want to write the mesh to disk since we use this data in
      // the multiscale analysis
      if (!print_field_frequency)
      {
        copyData(u_h, u_d);
        copyData(v_h, v_d);
        copyData(a_h, a_d);
        copyData(f_int_h, f_int_d);
        copyData(f_ext_d, f_ext_h);
      }
      if (print_history_frequency && (n_step % print_history_frequency))
      {
        writer->writeHistoryData(itr_prev + n_step, W_int, W_ext, W_damp, W_kin, current_time);
      }
      itr_prev = n_step;
      BIO_V3(std::cout << "Microscale has successfully completed in " << n_step
                       << " timesteps.\n";)
      return true;
    }
    void setDispBC(int nfixed, int * dof, double * init_values, double * values)
    {
      disp_nfixed = nfixed;
      disp_boundary_dof_h = createHostIntArrayFromRaw(dof, nfixed);
      disp_boundary_values_h = createHostDoubleArrayFromRaw(values, nfixed);
      disp_boundary_init_values_h =
          createHostDoubleArrayFromRaw(init_values, nfixed);
    }
    // this function must be marked KOKKOS_INLINE_FUNCTION because it needs to
    // be run on the device
    KOKKOS_INLINE_FUNCTION
    double getLinearReactionForce(double orig_length,
                                  double length,
                                  double elastic_modulus,
                                  double area)
    {
      // abaqus ...
      // double length_ratio = length / orig_length;
      // double log_strain = log(length_ratio);
      // return elastic_modulus * area * log_strain / length_ratio;
      double length_ratio = length / orig_length;
      double green_strain = 1.0 / 2.0 * (length_ratio * length_ratio - 1);
      return length_ratio * elastic_modulus * area * green_strain;
    }
    protected:
    const std::string & analysis_name;
    Amplitude * amp;
    apf::Mesh2 * mesh;
    apf::Field * u_field;
    apf::Field * v_field;
    apf::Field * a_field;
    apf::Field * f_int_field;
    apf::Field * f_ext_field;
    apf::Field * nodalMass;
    apf::Field * coordinate_field;
    bool mass_field_initialized;
    double fiber_density;
    double fiber_elastic_modulus;
    double fiber_area;
    // arrays that will be unmanaged memory as they come from apf
    // these all will have a postfix _arr
    double * u_arr;
    double * v_arr;
    double * a_arr;
    double * f_int_arr;
    double * f_ext_arr;
    double * mass_matrix_arr;
    int * connectivity_arr;
    double * coords_arr;
    // Now we need host and device copies of these arrays.
    H_RWDA u_h;
    H_RWDA v_h;
    H_RWDA a_h;
    H_RWDA mass_matrix_h;
    H_RWIA connectivity_h;
    H_RWDA coords_h;
    // now we need device versions of these arrays
    D_RWDA u_d;
    D_RWDA v_d;
    D_RWDA a_d;
    D_RWDA mass_matrix_d;
    D_RWIA connectivity_d;
    D_RWDA coords_d;
    double visc_damp_coeff;
    bool print_field_by_num_frames;
    unsigned long print_field_frequency;
    unsigned long print_history_frequency;
    double total_time;
    double crit_time_scale_factor;
    int disp_nfixed;
    H_RWIA disp_boundary_dof_h;
    H_RWDA disp_boundary_values_h;
    H_RWDA disp_boundary_init_values_h;
    D_RWIA disp_boundary_dof_d;
    D_RWDA disp_boundary_values_d;
    D_RWDA disp_boundary_init_values_d;
    double energy_check_eps;
    int nelem;
    int nnds;
    int ndof;
    /// current time step
    unsigned long n_step;
    double dt_crit;
    double current_time;
    double t_npone;
    double t_nphalf;
    double dt_nphalf;
    double W_int;
    double W_ext;
    double W_damp;
    double W_kin;
    D_RWDA du_d;
    D_RWDA f_int_d;
    H_RWDA f_int_h;
    D_RWDA f_int_last_d;
    D_RWDA f_ext_d;
    H_RWDA f_ext_h;
    D_RWDA f_ext_last_d;
    D_RWDA f_damp_d;
    D_RWDA f_damp_last_d;
    D_RWDA f_d;
    D_RWDA l0_d;
    D_RWDA l_d;
    D_RWDA current_coords_d;
    ExplicitOutputWriter * writer;
    bool own_writer;
    protected:
    template <typename D>
    void deleteArray(D data)
    {
      static_cast<T &>(*this).deleteArray_(data);
    }
    void zeroHostData(H_RWDA data)
    {
      static_cast<T &>(*this).zeroHostData_(data);
    }
    void zeroDeviceData(D_RWDA data)
    {
      static_cast<T &>(*this).zeroDeviceData_(data);
    }
    H_RWDA createHostDoubleArray(int size)
    {
      return static_cast<T &>(*this).createHostDoubleArray_(size);
    }
    H_RWIA createHostIntArray(int size)
    {
      return static_cast<T &>(*this).createHostIntArray_(size);
    }
    D_RWDA createDeviceDoubleArray(int size)
    {
      return static_cast<T &>(*this).createDeviceDoubleArray_(size);
    }
    D_RWIA createDeviceIntArray(int size)
    {
      return static_cast<T &>(*this).createDeviceIntArray_(size);
    }
    D_RWIA createDeviceIntMirrorArray(H_RWIA hostData)
    {
      return static_cast<T &>(*this).createDeviceIntMirrorArray_(hostData);
    }
    D_RWDA createDeviceDoubleMirrorArray(H_RWDA hostData)
    {
      return static_cast<T &>(*this).createDeviceDoubleMirrorArray_(hostData);
    }
    // this is the default which is a shallow copy
    // kokkos will specialize this to a deep copy
    template <typename TO, typename FROM>
    void copyData_(TO & toArray, FROM fromArray)
    {
      toArray = fromArray;
    }
    template <typename TO, typename FROM>
    void copyData(TO & toArray, FROM fromArray)
    {
      static_cast<T &>(*this).copyData_(toArray, fromArray);
    }
    // note
    // note the underlying host array will not assume ownership
    // of the data.
    H_RWDA createHostDoubleArrayFromRaw_(double * raw_array, int)
    {
      return raw_array;
    }
    H_RWDA createHostDoubleArrayFromRaw(double * raw_array, int size)
    {
      return static_cast<T &>(*this).createHostDoubleArrayFromRaw_(raw_array,
                                                                   size);
    }
    H_RWIA createHostIntArrayFromRaw_(int * raw_array, int)
    {
      return raw_array;
    }
    H_RWIA createHostIntArrayFromRaw(int * raw_array, int size)
    {
      return static_cast<T &>(*this).createHostIntArrayFromRaw_(raw_array,
                                                                size);
    }
    void computeMassMatrix(int /* unused */,
                           apf::Mesh * mesh,
                           apf::Field * nodalMass,
                           double density,
                           double fiber_area)
    {
      // to get analytical solution use third order integration
      bio::MassIntegrator massInt(nodalMass, density, fiber_area, 3,
                                  bio::MassLumpType::RowSum);
      massInt.process(mesh, 1);
    }
    void getCurrentCoords(int ndof,
                          D_RODA coords,
                          D_RODA u,
                          D_RWDA current_coords)
    {
      static_cast<T &>(*this).getCurrentCoords_(ndof, coords, u,
                                                current_coords);
    }
    void getElementLengths(int nelem,
                           D_RODA coords,
                           D_ROIA connectivity,
                           D_RWDA l0)
    {
      static_cast<T &>(*this).getElementLengths_(nelem, coords, connectivity,
                                                 l0);
    }
    double getForces(int nelem,
                     int ndof,
                     D_RODA l0,
                     D_RODA l,
                     D_RODA v,
                     D_RODA current_coords,
                     D_ROIA connectivity,
                     D_RODA mass_matrix,
                     double visc_damp_coeff,
                     D_RWDA f_int,
                     D_RWDA f_int_last,
                     D_RWDA f_ext,
                     D_RWDA f_ext_last,
                     D_RWDA f_damp,
                     D_RWDA f_damp_last,
                     D_RWDA f,
                     double & residual)
    {
      return static_cast<T &>(*this).getForces_(
          nelem, ndof, fiber_elastic_modulus, fiber_area, fiber_density, l0, l,
          v, current_coords, connectivity, mass_matrix, visc_damp_coeff, f_int,
          f_int_last, f_ext, f_ext_last, f_damp, f_damp_last, f, residual);
    }
    void updateAcceleration(int ndof, D_RODA mass_matrix, D_RODA f, D_RWDA a)
    {
      static_cast<T &>(*this).updateAcceleration_(ndof, mass_matrix, f, a);
    }
    void updateVelocity(int ndof, D_RODA a, double dt, D_RWDA v)
    {
      static_cast<T &>(*this).updateVelocity_(ndof, a, dt, v);
    }
    void updateAccelVel(int ndof, double dt, D_RODA mass_matrix, D_RODA f,
                        D_RWDA a, D_RWDA v)
    {
      static_cast<T &>(*this).updateAccelVel_(ndof, dt, mass_matrix, f, a, v);
    }
    void applyDispBC(int nfixed,
                     double prev_t,
                     double t,
                     D_ROIA dof,
                     D_RODA init_values,
                     D_RODA values,
                     const Amplitude * amp,
                     D_RWDA u,
                     D_RWDA du)
    {
      double amp_t = amp->operator()(t);
      double amp_prev_t = amp->operator()(prev_t);
      static_cast<T &>(*this).applyDispBC_(nfixed, amp_t, amp_prev_t, dof,
                                           init_values, values, u, du);
    }
    void applyVelBC(int nfixed,
                    double t,
                    D_ROIA dof,
                    D_RODA values,
                    const Amplitude * amp,
                    D_RWDA a)
    {
      double amp_t = amp->derivative(t);
      static_cast<T &>(*this).applyVelBC_(nfixed, amp_t, dof, values, a);
    }
    void applyAccelBC(int nfixed,
                      double t,
                      D_ROIA dof,
                      D_RODA values,
                      const Amplitude * amp,
                      D_RWDA v)
    {
      double amp_t = amp->secondDerivative(t);
      static_cast<T &>(*this).applyAccelBC_(nfixed, amp_t, dof, values, v);
    }
    // set the force on the boundary to be 0, and the external forces to be
    // opotite of the internal forces
    void fixBoundaryForces(int nfixed,
                           D_ROIA dof,
                           D_RODA mass_matrix,
                           D_RODA a,
                           D_RWDA f_int,
                           D_RWDA f_ext,
                           D_RWDA f_damp,
                           D_RWDA f)
    {
      static_cast<T &>(*this).fixBoundaryForces_(nfixed, dof, mass_matrix, a,
                                                 f_int, f_ext, f_damp, f);
    }
    void applyAccelVelBC(int nfixed, double t, double visc_damp_coeff, D_ROIA dof, D_RODA values, D_RODA mass_matrix, const Amplitude * amp,
                         D_RWDA a, D_RWDA v, D_RWDA f_int, D_RWDA f_ext, D_RWDA f_damp, D_RWDA f)
    {
      double v_amp = amp->derivative(t);
      double a_amp = amp->secondDerivative(t);
      static_cast<T&>(*this).applyAccelVelBC_(nfixed, a_amp, v_amp, visc_damp_coeff,
                                              dof, values, mass_matrix,
                                              a, v, f_int, f_ext, f_damp, f);
    }
    void updateDisplacement(int ndof, D_RODA v, double dt, D_RWDA u, D_RWDA du)
    {
      static_cast<T &>(*this).updateDisplacement_(ndof, v, dt, u, du);
    }
    void computeEnergies(int ndof,
                         D_RODA du,
                         D_RODA f_int_last,
                         D_RODA f_int,
                         D_RODA f_ext_last,
                         D_RODA f_ext,
                         D_RODA f_damp_last,
                         D_RODA f_damp,
                         double & W_int,
                         double & W_ext,
                         double & W_damp)
    {
      static_cast<T &>(*this).computeEnergies_(ndof, du, f_int_last, f_int,
                                               f_ext_last, f_ext, f_damp_last,
                                               f_damp, W_int, W_ext, W_damp);
    }
    double computeKineticEnergy(int ndof, D_RODA mass_matrix, D_RODA v)
    {
      return static_cast<T &>(*this).computeKineticEnergy_(ndof, mass_matrix,
                                                           v);
    }
    bool checkEnergyBalance_(double W_kin,
                             double W_int,
                             double W_ext,
                             double eps)
    {
      double energy_bal = W_kin + W_int + W_damp - W_ext;
      return (fabs(energy_bal) <=
              eps * std::max(std::max(fabs(W_ext),
                                      std::max(fabs(W_int), fabs(W_kin))),
                             fabs(W_damp)));
    }
    bool checkEnergyBalance(double W_kin,
                            double W_int,
                            double W_ext,
                            double eps)
    {
      return static_cast<T &>(*this).checkEnergyBalance_(W_kin, W_int, W_ext,
                                                         eps);
    }
    void fence() { static_cast<T &>(*this).fence_(); }
    protected:
    // this protects us from accidentally calling CRTP with the wrong
    // instantiation
    ExplicitAnalysisBase(apf::Mesh2 * mesh,
                         double total_time,
                         double fiber_elastic_modulus,
                         double fiber_area,
                         double fiber_density,
                         Amplitude * amp,
                         const std::string & analysis_name,
                         double visc_damp_coeff = 0.01,
                         unsigned long print_history_frequency = 100000,
                         unsigned long print_field_frequency = 1000000,
                         bool print_field_by_num_frames = false,
                         double crit_time_scale_factor = 0.8,
                         double energy_check_eps = 1E-2,
                         apf::Field * coordinate_field = NULL,
                         apf::Field * u_field = NULL,
                         apf::Field * v_field = NULL,
                         apf::Field * a_field = NULL,
                         apf::Field * f_int_field = NULL,
                         apf::Field * f_ext_field = NULL,
                         apf::Field * mass_field = NULL,
                         ExplicitOutputWriter * writer = NULL)
        : analysis_name(analysis_name)
        , amp(amp)
        , mesh(mesh)
        , u_field(u_field)
        , v_field(v_field)
        , a_field(a_field)
        , f_int_field(f_int_field)
        , f_ext_field(f_ext_field)
        , nodalMass(mass_field)
        , coordinate_field(coordinate_field)
        , mass_field_initialized(false)
        , fiber_density(fiber_density)
        , fiber_elastic_modulus(fiber_elastic_modulus)
        , fiber_area(fiber_area)
        , visc_damp_coeff(visc_damp_coeff)
        , print_field_by_num_frames(print_field_by_num_frames)
        , print_field_frequency(print_field_frequency)
        , print_history_frequency(print_history_frequency)
        , total_time(total_time)
        , crit_time_scale_factor(crit_time_scale_factor)
        , energy_check_eps(energy_check_eps)
        , writer(writer)
        , own_writer(false)
    {
      if (coordinate_field == NULL)
      {
        coordinate_field = mesh->getCoordinateField();
      }
      if (u_field == NULL)
      {
        u_field = apf::createLagrangeField(mesh, "u", apf::VECTOR, 1);
        apf::zeroField(u_field);
      }
      if (v_field == NULL)
      {
        v_field = apf::createLagrangeField(mesh, "v", apf::VECTOR, 1);
        apf::zeroField(v_field);
      }
      if (a_field == NULL)
      {
        a_field = apf::createLagrangeField(mesh, "a", apf::VECTOR, 1);
        apf::zeroField(a_field);
      }
      if (f_int_field == NULL)
      {
        f_int_field = apf::createLagrangeField(mesh, "f", apf::VECTOR, 1);
        apf::zeroField(f_int_field);
      }
      if (f_ext_field == NULL)
      {
        f_ext_field = apf::createLagrangeField(mesh, "f_ext", apf::VECTOR, 1);
        apf::zeroField(f_ext_field);
      }
      if (mass_field == NULL)
      {
        mass_field_initialized = false;
        nodalMass = apf::createLagrangeField(mesh, "nodalMass", apf::SCALAR, 1);
        apf::zeroField(nodalMass);
      }
      if (writer == NULL)
      {
        writer = new ExplicitOutputWriter(mesh,
                 (analysis_name).c_str(),
                 (analysis_name + ".pvd").c_str());
          own_writer = true;
      }
      // convert apf fields to arrays
      apf::freeze(u_field);
      apf::freeze(v_field);
      apf::freeze(a_field);
      apf::freeze(f_int_field);
      apf::freeze(f_ext_field);
      apf::freeze(nodalMass);
      u_arr = apf::getArrayData(u_field);
      v_arr = apf::getArrayData(v_field);
      a_arr = apf::getArrayData(a_field);
      f_int_arr = apf::getArrayData(f_int_field);
      f_ext_arr = apf::getArrayData(f_ext_field);
      mass_matrix_arr = apf::getArrayData(nodalMass);
      int etype;
      apf::destruct(mesh, connectivity_arr, nelem, etype, 1);
      assert(etype == apf::Mesh::EDGE);
      // apf::extractCoords(mesh, coords_arr, nnds);
      // TODO Check to see if we can accomplish the same thing by
      // just freezing the xpy field. My concern is that that field
      // is a function field and may behave badly...
      extractCoordinateArray(mesh, coordinate_field, coords_arr, nnds);
      ndof = nnds * 3;
      n_step = 0;
      dt_crit = std::numeric_limits<double>::max();
      current_time = 0;
      t_npone = 0;
      t_nphalf = 0;
      dt_nphalf = 0;
      W_int = 0;
      W_ext = 0;
      W_damp = 0;
      W_kin = 0;
      disp_nfixed = 0;
      // create the host arrays from the raw pointer data
      u_h = createHostDoubleArrayFromRaw(u_arr, ndof);
      v_h = createHostDoubleArrayFromRaw(v_arr, ndof);
      a_h = createHostDoubleArrayFromRaw(a_arr, ndof);
      f_int_h = createHostDoubleArrayFromRaw(f_int_arr, ndof);
      f_ext_h = createHostDoubleArrayFromRaw(f_ext_arr, ndof);
      mass_matrix_h = createHostDoubleArrayFromRaw(mass_matrix_arr, nnds);
      connectivity_h = createHostIntArrayFromRaw(connectivity_arr, 2 * nelem);
      coords_h = createHostDoubleArrayFromRaw(coords_arr, ndof);
      u_d = createDeviceDoubleMirrorArray(u_h);
      v_d = createDeviceDoubleMirrorArray(v_h);
      a_d = createDeviceDoubleMirrorArray(a_h);
      f_int_d = createDeviceDoubleMirrorArray(f_int_h);
      f_ext_d = createDeviceDoubleMirrorArray(f_int_h);
      mass_matrix_d = createDeviceDoubleMirrorArray(mass_matrix_h);
      connectivity_d = createDeviceIntMirrorArray(connectivity_h);
      coords_d = createDeviceDoubleMirrorArray(coords_h);
    }
    ~ExplicitAnalysisBase() {
      if(own_writer)
       delete writer; 
      delete [] coords_arr;
      delete [] connectivity_arr;
    }
    friend T;
  };
  using SH_RODA = const double *;
  using SH_RWDA = double *;
  using SH_ROIA = const int *;
  using SH_RWIA = int *;
  using SD_RODA = const double *;
  using SD_RWDA = double *;
  using SD_ROIA = const int *;
  using SD_RWIA = int *;
  // note that for the serial case, all of the host and device, and read only
  // arrays have the same types. The need for all of these types as part
  // of the class definition is a bit clunky, but it's better than needing
  // to provide types for each function indivitually.
  class ExplicitAnalysisSerial
      : public ExplicitAnalysisBase<ExplicitAnalysisSerial,
                                    SH_RODA,
                                    SH_RWDA,
                                    SH_ROIA,
                                    SH_RWIA,
                                    SD_RODA,
                                    SD_RWDA,
                                    SD_ROIA,
                                    SD_RWIA,
                                    SH_RODA,
                                    SH_RODA>
  {
    public:
    ExplicitAnalysisSerial(apf::Mesh2 * mesh,
                           double total_time,
                           double fiber_elastic_modulus,
                           double fiber_area,
                           double fiber_density,
                           Amplitude * amp,
                           const std::string & analysis_name,
                           double visc_damp_coeff = 0.01,
                           unsigned long print_history_frequency = 100000,
                           unsigned long print_field_frequency = 1000000,
                           bool print_field_by_num_frames = false,
                           double crit_time_scale_factor = 0.8,
                           double energy_check_eps = 1E-2,
                           apf::Field * coordinate_field = NULL,
                           apf::Field * u_field = NULL,
                           apf::Field * v_field = NULL,
                           apf::Field * a_field = NULL,
                           apf::Field * f_int_field = NULL,
                           apf::Field * f_ext_field = NULL,
                           apf::Field * mass_field = NULL,
                           ExplicitOutputWriter * writer = NULL)
        : ExplicitAnalysisBase(mesh,
                               total_time,
                               fiber_elastic_modulus,
                               fiber_area,
                               fiber_density,
                               amp,
                               analysis_name,
                               visc_damp_coeff,
                               print_history_frequency,
                               print_field_frequency,
                               print_field_by_num_frames,
                               crit_time_scale_factor,
                               energy_check_eps,
                               coordinate_field,
                               u_field,
                               v_field,
                               a_field,
                               f_int_field,
                               f_ext_field,
                               mass_field,
                               writer)
    {
      du_d = createDeviceDoubleArray(ndof);
      f_d = createDeviceDoubleArray(ndof);
      f_int_last_d = createDeviceDoubleArray(ndof);
      // f_ext_d = createDeviceDoubleArray(ndof);
      f_ext_last_d = createDeviceDoubleArray(ndof);
      f_damp_d = createDeviceDoubleArray(ndof);
      f_damp_last_d = createDeviceDoubleArray(ndof);
      l0_d = createDeviceDoubleArray(nelem);
      l_d = createDeviceDoubleArray(nelem);
      current_coords_d = createDeviceDoubleArray(ndof);
      zeroDeviceData(f_int_last_d);
      // zeroDeviceData(f_ext_d, ndof);
      zeroDeviceData(f_ext_last_d);
      zeroDeviceData(f_damp_d);
      zeroDeviceData(f_damp_last_d);
      zeroDeviceData(f_d);
    }
    ~ExplicitAnalysisSerial()
    {
      deleteArray(du_d);
      deleteArray(f_int_last_d);
      // deleteArray(f_ext_d);
      deleteArray(f_ext_last_d);
      deleteArray(f_damp_d);
      deleteArray(f_damp_last_d);
      deleteArray(f_d);
      deleteArray(l0_d);
      deleteArray(l_d);
      deleteArray(current_coords_d);
    }

    public:
    void zeroDeviceData_(double * data) { std::fill_n(data, ndof, 0); }
    void zeroHostData_(double * data) { std::fill_n(data, ndof, 0); }
    double * createHostDoubleArray_(int size) { return new double[size]; }
    int * createHostIntArray_(int size) { return new int[size]; }
    double * createDeviceDoubleArray_(int size) { return new double[size]; }
    int * createDeviceIntArray_(int size) { return new int[size]; }
    template <typename D>
    void deleteArray_(D data)
    {
      delete[] data;
    }
    SD_RWDA createDeviceDoubleMirrorArray_(SH_RWDA hostData)
    {
      return hostData;
    }
    SD_RWIA createDeviceIntMirrorArray_(SH_RWIA hostData) { return hostData; }
    void getCurrentCoords_(int ndof,
                           const double * coords,
                           const double * u,
                           double * current_coords)
    {
      for (int i = 0; i < ndof; ++i)
      {
        current_coords[i] = coords[i] + u[i];
      }
    }
    void getElementLengths_(int nelem,
                            const double * coords,
                            const int * connectivity,
                            double * l0)
    {
      int n1, n2;
      double x1, x2, x3;
      for (int i = 0; i < nelem; ++i)
      {
        n1 = connectivity[i * 2];
        n2 = connectivity[i * 2 + 1];
        x1 = coords[n2 * 3] - coords[n1 * 3];
        x2 = coords[n2 * 3 + 1] - coords[n1 * 3 + 1];
        x3 = coords[n2 * 3 + 2] - coords[n1 * 3 + 2];
        l0[i] = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
      }
    }
    double getForces_(int nelem,
                      int ndof,
                      double fiber_elastic_modulus,
                      double fiber_area,
                      double fiber_density,
                      const double * l0,
                      const double * l,
                      const double * v,
                      const double * current_coords,
                      const int * connectivity,
                      const double * mass_matrix,
                      double visc_damp_coeff,
                      double * f_int,
                      double * f_int_last,
                      double * f_ext,
                      double * f_ext_last,
                      double * f_damp,
                      double * f_damp_last,
                      double * f,
                      double & residual)
    {
      double dt = std::numeric_limits<double>::max();
      double frc;
      // element normal vectors
      double elem_nrm_1, elem_nrm_2, elem_nrm_3;
      int n1, n2;
      double sound_speed = sqrt(fiber_elastic_modulus / fiber_density);
      assert(sound_speed > 1E-15);
      double dt_crit_elem;
      // swap force arrays and zero internal forces
      for (int i = 0; i < ndof; ++i)
      {
        f_int_last[i] = f_int[i];
        f_int[i] = 0;
        f_ext_last[i] = f_ext[i];
        // f_ext[i] = 0;
        f_damp_last[i] = f_damp[i];
        f_damp[i] = visc_damp_coeff * mass_matrix[i / 3] * v[i];
      }
      // set the internal forces
      for (int i = 0; i < nelem; ++i)
      {
        assert(l0[i] > 0);
        frc = getLinearReactionForce(l0[i], l[i], fiber_elastic_modulus,
                                     fiber_area);
        n1 = connectivity[i * 2];
        n2 = connectivity[i * 2 + 1];
        assert(l[i] > 0);
        elem_nrm_1 = (current_coords[n2 * 3] - current_coords[n1 * 3]) / l[i];
        elem_nrm_2 =
            (current_coords[n2 * 3 + 1] - current_coords[n1 * 3 + 1]) / l[i];
        elem_nrm_3 =
            (current_coords[n2 * 3 + 2] - current_coords[n1 * 3 + 2]) / l[i];
        f_int[n1 * 3] -= frc * elem_nrm_1;
        f_int[n1 * 3 + 1] -= frc * elem_nrm_2;
        f_int[n1 * 3 + 2] -= frc * elem_nrm_3;
        f_int[n2 * 3] += frc * elem_nrm_1;
        f_int[n2 * 3 + 1] += frc * elem_nrm_2;
        f_int[n2 * 3 + 2] += frc * elem_nrm_3;
        // possibly swap these rather than copying data
        dt_crit_elem = l[i] / sound_speed;
        if (dt_crit_elem < dt)
        {
          dt = dt_crit_elem;
        }
      }
      residual = 0;
      for (int i = 0; i < ndof; ++i)
      {
        double local_residual = f_ext[i]-f_int[i];
        residual +=  local_residual*local_residual;
        //f[i] = f_ext[i] - (f_int[i] + f_damp[i]);
        f[i] = local_residual - f_damp[i];
      }
      residual = sqrt(residual);
      return dt;
    }
    void updateAcceleration_(int ndof,
                             const double * mass_matrix,
                             const double * f,
                             double * a)
    {
      for (int i = 0; i < ndof; ++i)
      {
        a[i] = (1.0 / mass_matrix[i / 3]) * f[i];
      }
    }
    void updateVelocity_(int ndof, const double * a, double dt, double * v)
    {
      for (int i = 0; i < ndof; ++i)
      {
        v[i] = v[i] + dt * a[i];
      }
    }
    void updateAccelVel_(int ndof, double dt, const double * mass_matrix, const double* f, 
                         double * a, double * v)
    {
      // 2*ndof loads, 2*ndof writes
      // 2*ndof multiply, 2*ndof divide
      for(int i=0; i<ndof; ++i)
      {
        double a_local = (1.0 / mass_matrix[i / 3]) * f[i];  
        a[i] = a_local;
        v[i] += dt * a_local;
      }
    }
    void updateDisplacement_(int ndof,
                             const double * v,
                             double dt,
                             double * u,
                             double * du)
    {
      for (int i = 0; i < ndof; ++i)
      {
        double du_local = dt * v[i];
        du[i] = du_local;
        u[i] += du_local;
      }
    }
    void applyDispBC_(int nfixed,
                      double amp_t,
                      double amp_prev_t,
                      const int * dof,
                      const double * init_values,
                      const double * values,
                      double * u,
                      double * du)
    {
      for (int i = 0; i < nfixed; ++i)
      {
        u[dof[i]] = values[i] * amp_t + init_values[i];
        du[dof[i]] = values[i] * (amp_t - amp_prev_t);
        // assert(fabs(u[dof[i]]) < fabs(values[i]) ||
        // isClose(u[dof[i]],values[i]));
      }
    }
    void applyVelBC_(int nfixed,
                     double amp_t,
                     SD_ROIA dof,
                     SD_RODA values,
                     SD_RWDA v)
    {
      for (int i = 0; i < nfixed; ++i)
      {
        v[dof[i]] = values[i] * amp_t;
      }
    }
    void applyAccelBC_(int nfixed,
                       double amp_t,
                       SD_ROIA dof,
                       SD_RODA values,
                       SD_RWDA a)
    {
      for (int i = 0; i < nfixed; ++i)
      {
        a[dof[i]] = values[i] * amp_t;
      }
    }
    void fixBoundaryForces_(int nfixed,
                            const int * dof,
                            const double * mass_matrix,
                            const double * a,
                            const double * f_int,
                            double * f_ext,
                            double * f_damp,
                            double * f)
    {
      for (int i = 0; i < nfixed; ++i)
      {
        double f_inertial = mass_matrix[dof[i] / 3] * a[dof[i]];
        f_ext[dof[i]] = f_inertial + f_int[dof[i]] + f_damp[dof[i]];
        f[dof[i]] = f_inertial;
      }
    }
    void applyAccelVelBC_(int nfixed, double a_amp, double v_amp, double visc_damp_coeff,
                          SD_ROIA dof, SD_RODA values, SD_RODA mass_matrix,
                          SD_RWDA a, SD_RWDA v, SD_RWDA f_int, SD_RWDA f_ext, SD_RWDA f_damp,
                          SD_RWDA f)
    {
      //4*nfixed reads, 5*nfixed writes
      //5*nfixed multiply, 2*nfixed adds, 1 divide
      for(int i=0; i<nfixed; ++i)
      {
        int local_dof = dof[i];
        double value = values[i];
        double a_local = value * a_amp;
        double v_local = value * v_amp;
        double mass = mass_matrix[local_dof / 3];
        double f_damp_local = visc_damp_coeff * mass * v_local;
        double f_inertial = mass * a_local;
        a[local_dof] =  a_local;
        v[local_dof] = v_local;
        f_damp[local_dof] = f_damp_local;
        f_ext[local_dof] =
            f_inertial + f_int[local_dof] + f_damp_local;
        f[local_dof] = f_inertial;
      }
    }
    void computeEnergies_(int ndof,
                          SD_RODA du,
                          SD_RODA f_int_last,
                          SD_RODA f_int,
                          SD_RODA f_ext_last,
                          SD_RODA f_ext,
                          SD_RODA f_damp_last,
                          SD_RODA f_damp,
                          double & W_int,
                          double & W_ext,
                          double & W_damp)
    {
      double du_local;
      double W_int_local = 0;
      double W_ext_local = 0;
      double W_damp_local = 0;
      for (int i = 0; i < ndof; ++i)
      {
        du_local = du[i];
        W_int_local += du_local * (f_int_last[i] + f_int[i]);
        W_ext_local += du_local * (f_ext_last[i] + f_ext[i]);
        W_damp_local += du_local * (f_damp_last[i] + f_damp[i]);
      }
      W_int = W_int + 0.5 * W_int_local;
      W_ext = W_ext + 0.5 * W_ext_local;
      W_damp = W_damp + 0.5 * W_damp_local;
    }
    double computeKineticEnergy_(int ndof,
                                 const double * mass_matrix,
                                 const double * v)
    {
      double W_kin = 0;
      for (int i = 0; i < ndof; ++i)
      {
        W_kin += mass_matrix[i / 3] * v[i] * v[i];
      }
      return 0.5 * W_kin;
    }
    void fence_() {}
  };
  using exe_space = Kokkos::DefaultExecutionSpace;
  using KKD_RODA = Kokkos::View<const double *, exe_space::device_type>;
  using KKD_RWDA = Kokkos::View<double *, exe_space::device_type>;
  using KKD_ROIA = Kokkos::View<const int *, exe_space::device_type>;
  using KKD_RWIA = Kokkos::View<int *, exe_space::device_type>;
  using KKH_RODA = KKD_RODA::HostMirror;
  using KKH_RWDA = KKD_RWDA::HostMirror;
  using KKH_ROIA = KKD_ROIA::HostMirror;
  using KKH_RWIA = KKD_RWIA::HostMirror;
  using KKD_RARODA = Kokkos::View<const double *,
                                  exe_space::device_type,
                                  Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
  using KKD_RAROIA = Kokkos::View<const int *,
                                  exe_space::device_type,
                                  Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
  class ExplicitAnalysisKokkos
      : public ExplicitAnalysisBase<ExplicitAnalysisKokkos,
                                    KKH_RODA,
                                    KKH_RWDA,
                                    KKH_ROIA,
                                    KKH_RWIA,
                                    KKD_RODA,
                                    KKD_RWDA,
                                    KKD_ROIA,
                                    KKD_RWIA,
                                    KKD_RARODA,
                                    KKD_RARODA>
  {
    public:
    template <typename D>
    void deleteArray_(D /* unused */ )
    {
    }
    KKH_RWDA createHostDoubleArrayFromRaw_(double * raw_array, int size)
    {
      return Kokkos::View<double *, Kokkos::DefaultHostExecutionSpace,
                          Kokkos::MemoryUnmanaged>(raw_array, size);
    }
    KKH_RWIA createHostIntArrayFromRaw_(int * raw_array, int size)
    {
      return Kokkos::View<int *, Kokkos::DefaultHostExecutionSpace,
                          Kokkos::MemoryUnmanaged>(raw_array, size);
    }
    KKD_RWDA createDeviceDoubleArray_(int size)
    {
      return KKD_RWDA("dbl_arr", size);
    }
    KKD_RWIA createDeviceIntArray_(int size)
    {
      return KKD_RWIA("int_arr", size);
    }
    KKH_RWDA createHostDoubleArray_(int size)
    {
      return KKH_RWDA("dbl_arr", size);
    }
    KKH_RWIA createHostIntArray_(int size) { return KKH_RWIA("int_arr", size); }
    KKD_RWDA createDeviceDoubleMirrorArray_(KKH_RWDA hostData)
    {
      return Kokkos::create_mirror_view(exe_space(), hostData);
    }
    KKD_RWIA createDeviceIntMirrorArray_(KKH_RWIA hostData)
    {
      return Kokkos::create_mirror_view(exe_space(), hostData);
    }
    void zeroDeviceData(KKD_RWDA data, int n)
    {
      Kokkos::parallel_for("zeroDeviceData", n,
                           KOKKOS_LAMBDA(std::size_t i) { data(i) = 0; });
    }
    void getCurrentCoords_(int ndof,
                           KKD_RODA coords,
                           KKD_RODA u,
                           KKD_RWDA current_coords)
    {
      Kokkos::parallel_for("getCurrentCoords", ndof,
                           KOKKOS_LAMBDA(std::size_t i) {
                             current_coords(i) = coords(i) + u(i);
                           });
    }
    void getElementLengths_(int nelem,
                            KKD_RODA coords,
                            KKD_ROIA connectivity,
                            KKD_RWDA l0)
    {
      Kokkos::parallel_for(
          "getElementLengths", nelem, KOKKOS_LAMBDA(std::size_t i) {
            int n1 = connectivity(i * 2);
            int n2 = connectivity(i * 2 + 1);
            double x1 = coords(n2 * 3) - coords(n1 * 3);
            double x2 = coords(n2 * 3 + 1) - coords(n1 * 3 + 1);
            double x3 = coords(n2 * 3 + 2) - coords(n1 * 3 + 2);
            l0(i) = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
          });
    }
    double getForces_(int nelem,
                      int ndof,
                      double fiber_elastic_modulus,
                      double fiber_area,
                      double fiber_density,
                      KKD_RODA l0,
                      KKD_RODA l,
                      KKD_RODA v,
                      KKD_RODA current_coords,
                      KKD_ROIA connectivity,
                      KKD_RODA mass_matrix,
                      double visc_damp_coeff,
                      KKD_RWDA f_int,
                      KKD_RWDA f_int_last,
                      KKD_RWDA f_ext,
                      KKD_RWDA f_ext_last,
                      KKD_RWDA f_damp,
                      KKD_RWDA f_damp_last,
                      KKD_RWDA f,
                      double & residual)
    {
      double dt = std::numeric_limits<double>::max();
      // element normal vectors
      double sound_speed = sqrt(fiber_elastic_modulus / fiber_density);
      residual = 0;
      // swap force arrays and zero internal forces
      Kokkos::parallel_for(
          "getFoces--Loop1", ndof, KOKKOS_LAMBDA(std::size_t i) {
            f_int_last(i) = f_int(i);
            f_int(i) = 0;
            f_ext_last(i) = f_ext(i);
            // f_ext(i) = 0;
            f_damp_last(i) = f_damp(i);
            f_damp(i) = visc_damp_coeff * mass_matrix(i / 3) * v(i);
          });
      // set the internal forces
      // TODO test using a reducer functional rather than atomic adds,
      // also try doing this with hierarchical parallelism and loading
      // loading forces vector into the shared memory.
      Kokkos::Min<double> min_reducer(dt);
      Kokkos::parallel_reduce(
          "getForces--mainLoop", Kokkos::RangePolicy<exe_space>(0, nelem),
          KOKKOS_LAMBDA(std::size_t i, double & dt_crit_elem) {
            double local_l = l(i);
            // FIXME the force reactions are all messed up! need to figure out
            // how to get the struct data in here?
            double frc = getLinearReactionForce(
                l0(i), local_l, fiber_elastic_modulus, fiber_area);
            int n1 = connectivity(i * 2);
            int n2 = connectivity(i * 2 + 1);
            double elem_nrm_1 =
                (current_coords(n2 * 3) - current_coords(n1 * 3)) / local_l;
            double elem_nrm_2 =
                (current_coords(n2 * 3 + 1) - current_coords(n1 * 3 + 1)) /
                local_l;
            double elem_nrm_3 =
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
      Kokkos::parallel_reduce("getForces-Loop3",
                              Kokkos::RangePolicy<exe_space>(0,ndof),
                              KOKKOS_LAMBDA(std::size_t i, double & residual_update) {
        double local_residual = f_ext(i)-f_int(i);
        residual_update += local_residual*local_residual;
        f(i) = local_residual - f_damp(i);
      }, residual);
      residual = sqrt(residual);
      return dt;
    }
    void updateAcceleration_(int ndof,
                             KKD_RODA mass_matrix,
                             KKD_RODA f,
                             KKD_RWDA a)
    {
      Kokkos::parallel_for("updateAcceleration", ndof,
                           KOKKOS_LAMBDA(std::size_t i) {
                             a(i) = (1.0 / mass_matrix(i / 3)) * f(i);
                           });
    }
    void updateVelocity_(int ndof, KKD_RODA a, double dt, KKD_RWDA v)
    {
      Kokkos::parallel_for("updateVelocity", ndof,
                           KOKKOS_LAMBDA(std::size_t i) { v(i) += dt * a(i); });
    }
    void updateAccelVel_(int ndof, double dt, KKD_RODA mass_matrix, KKD_RODA f, KKD_RWDA a, KKD_RWDA v)
    {
      // 2*ndof loads, 2*ndof writes
      // 2*ndof multiply, 2*ndof divide
      Kokkos::parallel_for("updateAccelVel", ndof,
                           KOKKOS_LAMBDA(std::size_t i) {
                             double a_local = (1.0 / mass_matrix(i / 3)) * f(i);  
                             a(i) = a_local;
                             v(i) += dt * a_local;
                           });
    }
    void updateDisplacement_(int ndof,
                             KKD_RODA v,
                             double dt,
                             KKD_RWDA u,
                             KKD_RWDA du)
    {
      Kokkos::parallel_for("updateDisplacement", ndof,
                           KOKKOS_LAMBDA(std::size_t i) {
                             double du_local = dt * v(i);
                             du(i) = du_local;
                             u(i) += du_local;
                           });
    }
    void applyDispBC_(int nfixed,
                      double amp_t,
                      double amp_prev_t,
                      KKD_ROIA dof,
                      KKD_RODA init_values,
                      KKD_RODA values,
                      KKD_RWDA u,
                      KKD_RWDA du)
    {
      // FIXME here we have a pointer to a class which will fail in cuda!
      Kokkos::parallel_for("applyDispBC", nfixed, KOKKOS_LAMBDA(std::size_t i) {
        int dof_local = dof(i);
        double val = values(i);
        double du_local = val * amp_t;
        u(dof_local) = du_local + init_values(i);
        du(dof_local) = du_local - val * amp_prev_t;
      });
    }
    void applyVelBC_(int nfixed,
                     double amp_t,
                     KKD_ROIA dof,
                     KKD_RODA values,
                     KKD_RWDA v)
    {
      Kokkos::parallel_for("applyVelBC", nfixed, KOKKOS_LAMBDA(std::size_t i) {
        v(dof(i)) = values(i) * amp_t;
      });
    }
    void applyAccelBC_(int nfixed,
                       double amp_t,
                       KKD_ROIA dof,
                       KKD_RODA values,
                       KKD_RWDA a)
    {
      Kokkos::parallel_for("applyAccelBC", nfixed, KOKKOS_LAMBDA(std::size_t i) {
        a(dof(i)) = values(i) * amp_t;
      });
    }
    void fixBoundaryForces_(int nfixed,
                            KKD_ROIA dof,
                            KKD_RODA mass_matrix,
                            KKD_RODA a,
                            KKD_RWDA f_int,
                            KKD_RWDA f_ext,
                            KKD_RWDA f_damp,
                            KKD_RWDA f)
    {
      Kokkos::parallel_for(
          "fixBoundaryForce", nfixed, KOKKOS_LAMBDA(std::size_t i) {
            // if we apply a velocity, or displacement boundary
            // condition, the external force is the sum of the
            // internal, and inertial forces
            int local_dof = dof(i);
            double f_inertial = mass_matrix(local_dof / 3) * a(local_dof);
            f_ext(local_dof) =
                f_inertial + f_int(local_dof) + f_damp(local_dof);
            f(local_dof) = f_inertial;
          });
    }
    void applyAccelVelBC_(int nfixed, double a_amp, double v_amp, double visc_damp_coeff,
                          KKD_ROIA dof, KKD_RODA values, KKD_RODA mass_matrix,
                          KKD_RWDA a, KKD_RWDA v, KKD_RWDA f_int, KKD_RWDA f_ext, KKD_RWDA f_damp,
                          KKD_RWDA f)
    {
      //4*nfixed reads, 5*nfixed writes
      //5*nfixed multiply, 2*nfixed adds, 1 divide
      Kokkos::parallel_for("applyAccelVelBC", nfixed, KOKKOS_LAMBDA(std::size_t i) {
        int local_dof = dof(i);
        double value = values(i);
        double a_local = value * a_amp;
        double v_local = value * v_amp;
        double mass = mass_matrix(local_dof / 3);
        double f_damp_local = visc_damp_coeff * mass * v_local;
        double f_inertial = mass * a_local;
        a(local_dof) =  a_local;
        v(local_dof) = v_local;
        f_damp(local_dof) = f_damp_local;
        f_ext(local_dof) =
            f_inertial + f_int(local_dof) + f_damp_local;
        f(local_dof) = f_inertial;
      });
    }
    // Struct for reduction of array of values
    struct EnergySum
    {
      typedef double value_type[];
      typedef Kokkos::View<double *>::size_type size_type;
      size_type value_count;
      EnergySum(KKD_RODA du,
                KKD_RODA f_int_last,
                KKD_RODA f_int,
                KKD_RODA f_ext_last,
                KKD_RODA f_ext,
                KKD_RODA f_damp_last,
                KKD_RODA f_damp)
          : value_count(3)
          , du_(du)
          , f_int_last_(f_int_last)
          , f_int_(f_int)
          , f_ext_last_(f_ext_last)
          , f_ext_(f_ext)
          , f_damp_last_(f_damp_last)
          , f_damp_(f_damp)
      {
      }
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type i, value_type energies) const
      {
        double du_local = du_(i);
        energies[0] += du_local * (f_int_last_(i) + f_int_(i));
        energies[1] += du_local * (f_ext_last_(i) + f_ext_(i));
        energies[2] += du_local * (f_damp_last_(i) + f_damp_(i));
      }
      KOKKOS_INLINE_FUNCTION
      void join(volatile value_type dst, volatile value_type src)
      {
        dst[0] += src[0];
        dst[1] += src[1];
        dst[2] += src[2];
      }
      KOKKOS_INLINE_FUNCTION
      void init(value_type energies) const
      {
        energies[0] = 0;
        energies[1] = 0;
        energies[2] = 0;
      }
      private:
      Kokkos::View<const double *> du_;
      Kokkos::View<const double *> f_int_last_;
      Kokkos::View<const double *> f_int_;
      Kokkos::View<const double *> f_ext_last_;
      Kokkos::View<const double *> f_ext_;
      Kokkos::View<const double *> f_damp_last_;
      Kokkos::View<const double *> f_damp_;
    };
    void computeEnergies_(int ndof,
                          KKD_RODA du,
                          KKD_RODA f_int_last,
                          KKD_RODA f_int,
                          KKD_RODA f_ext_last,
                          KKD_RODA f_ext,
                          KKD_RODA f_damp_last,
                          KKD_RODA f_damp,
                          double & W_int,
                          double & W_ext,
                          double & W_damp)
    {
      double energies[3];
      Kokkos::parallel_reduce(ndof,
                              EnergySum(du, f_int_last, f_int, f_ext_last,
                                        f_ext, f_damp_last, f_damp),
                              energies);
      W_int = W_int + 0.5 * energies[0];
      W_ext = W_ext + 0.5 * energies[1];
      W_damp = W_damp + 0.5 * energies[2];
    }
    double computeKineticEnergy_(int ndof, KKD_RODA mass_matrix, KKD_RODA v)
    {
      double W_kin = 0;
      Kokkos::parallel_reduce(
          "computeKineticEnergy", Kokkos::RangePolicy<exe_space>(0, ndof),
          KOKKOS_LAMBDA(std::size_t i, double & W_kin_update) {
            W_kin_update += mass_matrix(i / 3) * v(i) * v(i);
          },
          W_kin);
      return 0.5 * W_kin;
    }
    template <typename TO, typename FROM>
    void copyData_(TO & toArray, FROM fromArray)
    {
      Kokkos::deep_copy(toArray, fromArray);
    }
    void fence_() { Kokkos::fence(); }
    ExplicitAnalysisKokkos(apf::Mesh2 * mesh,
                           double total_time,
                           double fiber_elastic_modulus,
                           double fiber_area,
                           double fiber_density,
                           Amplitude * amp,
                           const std::string & analysis_name,
                           double visc_damp_coeff = 0.01,
                           unsigned long print_history_frequency = 100000,
                           unsigned long print_field_frequency = 1000000,
                           bool print_field_by_num_frames = false,
                           double crit_time_scale_factor = 0.8,
                           double energy_check_eps = 5E-2,
                           apf::Field * coordinate_field = NULL,
                           apf::Field * u_field = NULL,
                           apf::Field * v_field = NULL,
                           apf::Field * a_field = NULL,
                           apf::Field * f_int_field = NULL,
                           apf::Field * f_ext_field = NULL,
                           apf::Field * mass_field = NULL,
                           ExplicitOutputWriter * writer = NULL)
        : ExplicitAnalysisBase(mesh,
                               total_time,
                               fiber_elastic_modulus,
                               fiber_area,
                               fiber_density,
                               amp,
                               analysis_name,
                               visc_damp_coeff,
                               print_history_frequency,
                               print_field_frequency,
                               print_field_by_num_frames,
                               crit_time_scale_factor,
                               energy_check_eps,
                               coordinate_field,
                               u_field,
                               v_field,
                               a_field,
                               f_int_field,
                               f_ext_field,
                               mass_field,
                               writer)
    {
      du_d = createDeviceDoubleArray(ndof);
      f_d = createDeviceDoubleArray(ndof);
      f_int_last_d = createDeviceDoubleArray(ndof);
      // f_ext_d = createDeviceDoubleArray(ndof);
      f_ext_last_d = createDeviceDoubleArray(ndof);
      f_damp_d = createDeviceDoubleArray(ndof);
      f_damp_last_d = createDeviceDoubleArray(ndof);
      l0_d = createDeviceDoubleArray(nelem);
      l_d = createDeviceDoubleArray(nelem);
      current_coords_d = createDeviceDoubleArray(ndof);
      zeroDeviceData(f_d, ndof);
      zeroDeviceData(f_int_last_d, ndof);
      // zeroDeviceData(f_ext_d, ndof);
      zeroDeviceData(f_ext_last_d, ndof);
      zeroDeviceData(f_damp_d, ndof);
      zeroDeviceData(f_damp_last_d, ndof);
    }
    ~ExplicitAnalysisKokkos()
    {
      deleteArray(du_d);
      deleteArray(f_d);
      deleteArray(f_int_last_d);
      // deleteArray(f_ext_d);
      deleteArray(f_ext_last_d);
      deleteArray(f_damp_d);
      deleteArray(f_damp_last_d);
      deleteArray(l0_d);
      deleteArray(l_d);
      deleteArray(current_coords_d);
    }
  };
}  // namespace bio
