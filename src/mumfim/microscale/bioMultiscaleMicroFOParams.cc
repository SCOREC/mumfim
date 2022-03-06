#include "bioMultiscaleMicroFOParams.h"
namespace amsi
{
  template <typename T>
  struct static_init
  {
    static_init() : val(), init(false) {}
    T val;
    bool init;
  };
  template <>
  MPI_Datatype mpi_type<mumfim::micro_fo_header>()
  {
    static static_init<MPI_Datatype> tp;
    if(!tp.init)
    {
      MPI_Type_contiguous(mumfim::NUM_HEADER_FIELDS,MPI_INTEGER,&tp.val);
      MPI_Type_commit(&tp.val);
      tp.init = true;
    }
    return tp.val;
  }
  template <>
  MPI_Datatype mpi_type<mumfim::micro_fo_params>()
  {
    static static_init<MPI_Datatype> tp;
    if(!tp.init)
    {
      MPI_Type_contiguous(mumfim::NUM_PARAM_FIELDS,MPI_DOUBLE,&tp.val);
      MPI_Type_commit(&tp.val);
      tp.init = true;
    }
    return tp.val;
  }
  template <>
  MPI_Datatype mpi_type<mumfim::micro_fo_init_data>()
  {
    static static_init<MPI_Datatype> tp;
    if(!tp.init)
    {
      MPI_Type_contiguous(12,MPI_DOUBLE,&tp.val);
      MPI_Type_commit(&tp.val);
      tp.init = true;
    }
    return tp.val;
  }
  template <>
  MPI_Datatype mpi_type<mumfim::micro_fo_data>()
  {
    static static_init<MPI_Datatype> tp;
    if(!tp.init)
    {
      MPI_Type_contiguous(9,MPI_DOUBLE,&tp.val);
      MPI_Type_commit(&tp.val);
      tp.init = true;
    }
    return tp.val;
  }
  template <>
  MPI_Datatype mpi_type<mumfim::micro_fo_result>()
  {
    static static_init<MPI_Datatype> tp;
    if(!tp.init)
    {
      MPI_Type_contiguous(6 + 3 + 36,MPI_DOUBLE,&tp.val);
      MPI_Type_commit(&tp.val);
      tp.init = true;
    }
    return tp.val;
  }
  template <>
  MPI_Datatype mpi_type<mumfim::micro_fo_step_result>()
  {
    static static_init<MPI_Datatype> tp;
    if(!tp.init)
    {
      MPI_Type_contiguous(9+9,MPI_DOUBLE,&tp.val);
      MPI_Type_commit(&tp.val);
      tp.init = true;
    }
    return tp.val;
  }
  template <>
  MPI_Datatype mpi_type<mumfim::micro_fo_solver>()
  {
    static static_init<MPI_Datatype> tp;
    if(!tp.init)
    {
      MPI_Type_contiguous(mumfim::NUM_MICRO_SOLVER_FIELDS,MPI_DOUBLE,&tp.val);
      MPI_Type_commit(&tp.val);
      tp.init = true;
    }
    return tp.val;
  }
  template <>
  MPI_Datatype mpi_type<mumfim::micro_fo_int_solver>()
  {
    static static_init<MPI_Datatype> tp;
    if(!tp.init)
    {
      MPI_Type_contiguous(mumfim::NUM_MICRO_SOLVER_INT_FIELDS,MPI_INTEGER,&tp.val);
      MPI_Type_commit(&tp.val);
      tp.init = true;
    }
    return tp.val;
  }
}
