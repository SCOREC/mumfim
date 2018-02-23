#ifndef BIO_MICROFOMULTISCALE_TYPES_H_
#define BIO_MICROFOMULTISCALE_TYPES_H_
#include <amsiMPI.h>
namespace bio
{
  enum header_fields
  {
    RVE_TYPE = 0,
    FIELD_ORDER = 1,
    ELEMENT_TYPE = 2,
    GAUSS_ID = 3,
    FIBER_REACTION = 4,
    IS_ORIENTED = 5,
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
    ORIENTATION_ALIGN = 8,
    NUM_PARAM_FIELDS = 9
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
  struct MicroFODatatypes
  {
    MPI_Datatype hdr;
    MPI_Datatype prm;
    MPI_Datatype ini;
    MPI_Datatype dat;
    MPI_Datatype rst;
    MicroFODatatypes()
    {
      MPI_Type_contiguous(NUM_HEADER_FIELDS,MPI_INTEGER,&hdr);
      MPI_Type_contiguous(NUM_PARAM_FIELDS,MPI_DOUBLE,&prm);
      MPI_Type_contiguous(4*3,MPI_DOUBLE,&ini);
      MPI_Type_contiguous(9,MPI_DOUBLE,&dat);
      MPI_Type_contiguous(4*3*6+9+9,MPI_DOUBLE,&rst);
      MPI_Type_commit(&ini);
      MPI_Type_commit(&dat);
      MPI_Type_commit(&rst);
    }
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
}
#endif
