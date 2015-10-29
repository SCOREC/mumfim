#ifndef MICROFOMULTISCALE_TYPES_H_
#define MICROFOMULTISCALE_TYPES_H_

#include <mpi.h>

  namespace Biotissue {

    // can only support a single element type atm...
/*
    struct micro_fo_init_data_lin_tet
    {
      double init_data[4 * 3];
    };

    struct micro_fo_init_data_lin_pyr
    {
      double init_dat[5 * 3];
    }

    union micro_fo_init_data_union
    {
      micro_fo_init_data_lin_tet lin_tet;
      micro_fo_init_data_lin_pyr lin_pyr;
    };
*/
    struct micro_fo_init_data
    {
      int element_type;
      int gauss_pt_id;
      double init_data[4 * 3];
    };

    struct micro_fo_data
    {
       double data[4 * 6]; // 6 for each vertex, assuming tets atm
    };

    struct micro_fo_result
    {
      // 8 * 3 * 6 + 9 = 153 (magic number, 8 verts in 'largest' element type)
      double data[4 * 3 * 6 + 9 + 3]; // 6 sigma values for each vertex 3 q values, assuming tets atm
    };
    
    class MicroFOMultiscaleDataTypes{
    
      public:
      MPI_Datatype micro_fo_init_data_type;    
      MPI_Datatype micro_fo_data_type;    
      MPI_Datatype micro_fo_result_type;

      void MultiscaleDataTypesMPICommit()
      {
        int init_block_lengths[] = {1,1,4*3};
        MPI_Aint init_block_disp[] = {0,sizeof(int),2*sizeof(int)};
        MPI_Datatype init_block_types[] = {MPI_INTEGER, MPI_INTEGER, MPI_DOUBLE};

        MPI_Type_create_struct(3,
			    init_block_lengths,
			    init_block_disp,
			    init_block_types,
			    &micro_fo_init_data_type);
	
        MPI_Type_contiguous(4*6,MPI_DOUBLE,&micro_fo_data_type);
	
        MPI_Type_contiguous(4*3*6+9+3,MPI_DOUBLE,&micro_fo_result_type);

        MPI_Type_commit(&micro_fo_init_data_type);
        MPI_Type_commit(&micro_fo_data_type);
        MPI_Type_commit(&micro_fo_result_type);
      }

    };


  }


#endif
