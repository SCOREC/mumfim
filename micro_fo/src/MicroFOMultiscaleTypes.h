#ifndef MICROFOMULTISCALE_TYPES_H_
#define MICROFOMULTISCALE_TYPES_H_
#include <mpi.h>
namespace bio
{
  enum header_fields
  {
    RVE_TYPE = 0,
    ELEMENT_TYPE = 1,
    GAUSS_ID = 2,
    FIBER_REACTION = 3,
    IS_ORIENTED = 4,
    NUM_HEADER_FIELDS = 5
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
    ORIENTATION_ALIGN = 8,
    NUM_PARAM_FIELDS = 9
  };
  struct micro_fo_params
  {
    double data[NUM_PARAM_FIELDS];
  };
  struct micro_fo_init_data
  {
<<<<<<< HEAD
    double init_data[9];
  };
  struct micro_fo_data
  {
    double data[9]; // 6 for each vertex, assuming tets atm
=======
    double init_data[4 * 3];
  };
  struct micro_fo_data
  {
    double data[4 * 6]; // 6 for each vertex, assuming tets atm
>>>>>>> develop
  };
  struct micro_fo_result
  {
    // 8 * 3 * 6 + 9 = 153 (magic number, 8 verts in 'largest' element type)
    // +9 at end to store orientation tensor information.
    double data[4 * 3 * 6 + 9 + 9]; // 6 sigma values for each vertex 3 q values, assuming tets atm
  };
  class MicroFOMultiscaleDataTypes{
  public:
    MPI_Datatype micro_fo_header_data_type;
    MPI_Datatype micro_fo_parameter_data_type;
    MPI_Datatype micro_fo_init_data_type;
    MPI_Datatype micro_fo_data_type;
    MPI_Datatype micro_fo_result_type;
    void MultiscaleDataTypesMPICommit()
    {
      MPI_Type_contiguous(NUM_HEADER_FIELDS,MPI_INTEGER,&micro_fo_header_data_type);
      MPI_Type_contiguous(NUM_PARAM_FIELDS,MPI_DOUBLE,&micro_fo_parameter_data_type);
      MPI_Type_contiguous(4*3,MPI_DOUBLE,&micro_fo_init_data_type);
      MPI_Type_contiguous(4*6,MPI_DOUBLE,&micro_fo_data_type);
      MPI_Type_contiguous(4*3*6+9+9,MPI_DOUBLE,&micro_fo_result_type);
      MPI_Type_commit(&micro_fo_init_data_type);
      MPI_Type_commit(&micro_fo_data_type);
      MPI_Type_commit(&micro_fo_result_type);
    }
  };
}
#endif
