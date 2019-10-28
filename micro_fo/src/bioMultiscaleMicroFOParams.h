#ifndef BIO_MULTISCALE_MICRO_FO_PARAMS_H_
#define BIO_MULTISCALE_MICRO_FO_PARAMS_H_
#include <amsiMPI.h>
namespace bio
{
  enum class MicroscaleType
  {
    NONE,
    FIBER_ONLY,
    FIBER_MATRIX,
    ISOTROPIC_NEOHOOKEAN,
    MICROSCALE_TYPE_COUNT
  };
  // TODO move orientation tensor fields to own structs
  // fiber_reaction is obsolete, remove it
  // keep the 2d orientation as a separate param because we might want to compute
  // multiple 2d orientations in the future
  enum header_fields
  {
    RVE_TYPE,
    RVE_DIR_TYPE,
    FIELD_ORDER,
    ELEMENT_TYPE,
    GAUSS_ID,
    COMPUTE_ORIENTATION_3D,
    COMPUTE_ORIENTATION_2D,
    NUM_HEADER_FIELDS
  };
  struct micro_fo_header
  {
    int data[NUM_HEADER_FIELDS];
  };
  enum fiber_param_fields
  {
    FIBER_RADIUS,
    VOLUME_FRACTION,
    YOUNGS_MODULUS,
    NONLINEAR_PARAM,
    LINEAR_TRANSITION,
    ORIENTATION_AXIS_X,
    ORIENTATION_AXIS_Y,
    ORIENTATION_AXIS_Z,
    NUM_PARAM_FIELDS
  };
  enum micro_solver_fields
  {
    MICRO_SOLVER_TOL,
    MICRO_CONVERGENCE_TOL,
    PREV_ITER_FACTOR,
    LOAD_TIME,
    HOLD_TIME,
    VISCOUS_DAMPING_FACTOR,
    CRITICAL_TIME_SCALE_FACTOR,
    ENERGY_CHECK_EPSILON,
    NUM_MICRO_SOLVER_FIELDS
  };
  enum micro_solver_int_fields
  {
    MICRO_SOLVER_TYPE,
    MAX_MICRO_CUT_ATTEMPT,
    MICRO_ATTEMPT_CUT_FACTOR,
    MAX_MICRO_ITERS,
    DETECT_OSCILLATION_TYPE,
    AMPLITUDE_TYPE,
    PRINT_HISTORY_FREQUENCY,
    PRINT_FIELD_FREQUENCY,
    PRINT_FIELD_BY_NUM_FRAMES,
    SERIAL_GPU_CUTOFF,
    NUM_MICRO_SOLVER_INT_FIELDS
  };
  struct micro_fo_solver
  {
    double data[NUM_MICRO_SOLVER_FIELDS];
  };
  struct micro_fo_int_solver
  {
    int data[NUM_MICRO_SOLVER_INT_FIELDS];
  };
  struct micro_fo_params
  {
    double data[NUM_PARAM_FIELDS];
  };
  // coordinates of the owning tet
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
    // 6 components of symmetric stress +
    // 3 components of Q + 
    // 36 components of stiffness (stress and strain symmetric)
    double data[6 + 3 + 36]; 
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
    MPI_Datatype mpi_type<bio::micro_fo_solver>();
  template <>
    MPI_Datatype mpi_type<bio::micro_fo_int_solver>();
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
