#ifndef SPARSKIT_H_
#define SPARSKIT_H_

#include <cmath>
#include <vector>

  extern "C" {

    void ilut_(int *,
	       double a[],
	       int ja[],
	       int ia[],
	       int *,
	       double *,
	       double *,
	       int *,
	       int *,
	       int *,
	       double *,
	       int *,
	       int *);
    
    void lusol_(int *,
		double a[],
		double *,
		double *,
		int *,
		int *);

    void pgmres_(int *,
		 int *,
		 double *,
		 double [],
		 double *,
		 double *,
		 int *,
		 int *,
		 double [],
		 int *,
		 int *,
		 double *,
		 int *,
		 int *,
		 int *);

    void amux_(int *,
	       double x[],
	       double y[],
	       double a[],
	       int ja[],
	       int ia[]);
    
  }

namespace Biotissue
{
  class SparskitBuffers
  {
  protected:
    int heuristic_length;
    std::vector<int> int_work_array;
    std::vector<int> rows;
    std::vector<int> cols;
    
    std::vector<double> double_work_array;
    std::vector<double> matrix;
  public:
    SparskitBuffers(int num_dofs);

    void zero();

    int matrixLength() {return heuristic_length;}
    int * intWorkBuffer() {return &int_work_array[0];}
    int * rowsBuffer() {return &rows[0];}
    int * colsBuffer() {return &cols[0];}
    double * doubleWorkBuffer() {return &double_work_array[0];}
    double * matrixBuffer() {return &matrix[0];}
  };
  
}


#endif
