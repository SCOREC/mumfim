#ifndef _CSR_MATRIX_H
#define _CSR_MATRIX_H

#include "listman.h"
#include <string>
#include <ostream>

class CSR_Vector;
class CSR_Matrix {

public :

  enum CSR_Matrix_Storage_t {HARWELL_BOEING_STORAGE, COMPRESSED_ROW_STORAGE, UNKNOWN_STORAGE};

  CSR_Matrix ();
  CSR_Matrix (const CSR_Matrix& in);  
  CSR_Matrix (const int n);
  CSR_Matrix & operator=(const CSR_Matrix& in);


  ~CSR_Matrix( );

  void Load (const CSR_Matrix& in);
 
  void ReAllocate(const int n, const int nbz);
  void SetSize(const int& n);

  void ExecuteReordering();

  void ReorderArray( double * array );
  void InverseReorderArray( double * array);

  void   AddMatrix(int i, int j, double val);
  double GetMatrix(int i, int j);
  void   ZeroMatrix  ( void );

  int      ReturnNumberOfNonZero ( void ) ;
  double * ReturnValPointer ( void ) ;
  int    * ReturnRowPointer  ( void ) ;
  int    * ReturnColumIndexPointer( void );

  void   EndOfAssembly();

  void   PutInCompressedRow( void );
  void   PutInHarwellBoeing( void );

  int GetNbUnknown();

  inline int ReturnStorageType()
  {return storage;};

  inline void SetMatrixReorderingToFalse( void ) 
  { MatrixWasReordered = false; }

  void OutputMatrixMatlabFormat(std::ostream &); 
  void cmkreord(int n, double *a, int *ja, int *ia, double *a0,
		int *ja0, int *ia0, int *init, int * iperm, int * mask,
		int * maskval, int * nlev, int * riord, int * levels);
  void sort_col(int n,double * a, int * ja, int *ia, int * iw, double * rw);
private :

  //to store the matrix
  List_T  *a_, *ai_, *ptr_, *jptr_;
  int NbLines;
  bool AllocationDone;
  CSR_Matrix_Storage_t storage;
  //for the reordering
  bool MatrixWasReordered;
  int     *permr_, *permp_, *rpermr_;


  void sort2(unsigned long n, double arr[], int ai[] , int aj [] );
  void deblign ( int nz , int *ptr , int *jptr , int *ai);
  int *ivector(long nl, long nh);
  void free_ivector(int *v, long nl, long nh);
  int  cmpij(int ai,int aj,int bi,int bj);


  int maskdeg(int *ja, int *ia, int *nod, int *mask, int *maskval);
  int rversp(int n, int *riord);
  int perphn(int n, int *ja, int *ia, int *init, int *iperm,
	     int * mask, int * maskval, int * nlev, int * riord, int * levels);
  void exchange(int n, int *iriord, int *iperm);
 
  int dperm(int nrow, double *a, int *ja, int * ia, double * ao,
	    int * jao, int * iao, int * perm, int *qperm, int * job);

  int cperm(int nrow, double * a, int *ja, int * ia, double * ao,
	    int * jao, int * iao, int * perm, int * job);

  void rperm(int nrow,double * a,int * ja,int * ia, double * ao, 
	     int * jao, int * iao, int * perm, int * job);
 
  int bfs(int n, int *ja,int *ia, int * nfirst, int * iperm, int * mask, 
	  int *maskval, int * riord,int *levels, int *nlev);

  int sort_irv(int * itmp, double *rtmp, int * n);

  int add_lvst(int * istart, int * iend, int * nlev, int * riord,
	       int * ja, int * ia, int * mask, int * maskval);

  void Allocate(const int n);
  void Clear(void);

};

void SPARSKIT_PRINT_POSTSCRIPT_ ( int fortran_file_descr , 
				  CSR_Matrix & mat );

void SPARSKIT_LINEAR_SOLVER_ ( const std::string &renumber, 
			       const std::string &prec, 
			       const std::string &accl, 
			       CSR_Matrix & mat , 
			       CSR_Vector & rhs, 
			       CSR_Vector &sol);
#endif








