#ifndef BIO_MICROFOMULTISCALE_TYPES_H_
#define BIO_MICROFOMULTISCALE_TYPES_H_
#include <amsiMPI.h>
namespace bio
{
  // TODO move orientation tensor fields to own structs
  // fiber_reaction is obsolete, remove it
  // keep the 2d orientation as a separate param because we might want to compute
  // multiple 2d orientations in the future
  enum header_fields
  {
    RVE_TYPE = 0,
    FIELD_ORDER = 1,
    ELEMENT_TYPE = 2,
    GAUSS_ID = 3,
    COMPUTE_ORIENTATION_3D = 4,
    COMPUTE_ORIENTATION_2D = 5,
    NUM_HEADER_FIELDS = 6
  };
  struct micro_fo_header
  {
    int data[NUM_HEADER_FIELDS];
  };
  enum fiber_param_fields
  {
    FIBER_RADIUS = 0,
    VOLUME_FRACTION = 1,
    YOUNGS_MODULUS = 2,
    NONLINEAR_PARAM = 3,
    LINEAR_TRANSITION = 4,
    ORIENTATION_AXIS_X = 5,
    ORIENTATION_AXIS_Y = 6,
    ORIENTATION_AXIS_Z = 7,
    NUM_PARAM_FIELDS = 8
  };
  struct micro_fo_params
  {
    double data[NUM_PARAM_FIELDS];
  };
  struct micro_fo_init_data
  {
    double init_data[4 * 3];
  };
  struct micro_fo_data
  {
    double data[9]; // deformation gradient
  };
  struct micro_fo_result
  {
    // 8 * 3 * 6 + 9 = 153 (magic number, 8 verts in 'largest' element type)
    // +9 at end to store orientation tensor information.
    double data[4 * 3 * 6 + 9 + 9]; // 6 sigma values for each vertex 3 q values, assuming tets atm
  };
  // data communicated at each step
  struct micro_fo_step_result
  {
    // 9 for 3D orientation tensor, and 9 for the 2D orientation tensor
    double data[9+9];
  };
}
namespace amsi
{
  template <>
    MPI_Datatype mpi_type<bio::micro_fo_header>();
  template <>
    MPI_Datatype mpi_type<bio::micro_fo_params>();
  template <>
    MPI_Datatype mpi_type<bio::micro_fo_init_data>();
  template <>
    MPI_Datatype mpi_type<bio::micro_fo_data>();
  template <>
    MPI_Datatype mpi_type<bio::micro_fo_result>();
  template <>
    MPI_Datatype mpi_type<bio::micro_fo_step_result>();
}
#endif
