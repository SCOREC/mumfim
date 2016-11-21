#include "mpi.h"
#include <string>
#include "CSR_Matrix.h"
#include "CSR_Vector.h"
#include <memory>


extern "C"
{
  /// running multiple solvers through fortran interface
  void pgmres_ (int *n, int *krylov_size, double *rhs, double *sol, double *wk, double *stoptest, int *, int *,  
		double *a, int *ja, int *ia,
		double *au, int *jau, int *ju, int *err);
  void runrc_ (int *n, double *rhs, double *sol, int *, double *, double *wk,   
	       double *guess, double *a, int *ja, int *ia,
	       double *au, int *jau, int *ju, int *err);
  void myplot_ (int *n, int *ja, int *ia, int *x, int *y);
  void ilu0_ (int *n, double *a, int *ja, int *ia, double *alu, int *jalu, int *ialu, int *x, int *y);
  void milu0_ (int *n, double *a, int *ja, int *ia, double *alu, int *jalu, int *ialu, int *x, int *y);
  void ilut_ (int *n, double *a, int *ja, int *ia, int *nbfill, double *droptol,
	      double *alu, int *jalu, int *ialu, int *nbElmILU, double *dw, int *iw, int *err);
  void lusol_(int *, double a[], double *, double *, int *, int *); // --LJZ 10/16/10 add back solve interface
}

void SPARSKIT_PRINT_POSTSCRIPT_ ( int fortran_file_descr , 
				  CSR_Matrix & mat )
{
  int nb_lines_matrix = mat.GetNbUnknown();
  int zero = 3;
  myplot_ (&nb_lines_matrix, 	 
	   mat.ReturnColumIndexPointer(), // matrix
	   mat.ReturnRowPointer(), // matrix
	   &fortran_file_descr, 
	   &zero); 
}

void SPARSKIT_LINEAR_SOLVER_ ( const std::string &renumber,
			       const std::string &prec, 
			       const std::string &accl, 			       
			       CSR_Matrix & mat , 
			       CSR_Vector & b, 
			       CSR_Vector & x)
{
  int ipar[256];
  double fpar[256];
  mat.PutInCompressedRow();


  int ierr = 0,i,j,k;
  double *  a    =  mat.ReturnValPointer();
  int *  ai    = mat.ReturnColumIndexPointer();
  int *  jptr      = mat.ReturnRowPointer();

     

  double *alu = 0;
  int    *jlu = 0;
  int    *ju = 0;

  double *a_rcmk    = 0;
  int *jptr_rcmk    = 0;
  int *ai_rcmk      = 0;


  int N = mat.GetNbUnknown();
  int nbz = mat.ReturnNumberOfNonZero () ;

  /* output the LHS matrix and RHS vector 
  FILE* LHS=fopen("LHSMatrix","w");
  fprintf(LHS, "%d\n", nbz);
  for(int i=0; i<nbz; i++)
      fprintf(LHS, "%0.16e\n", a[i]);
  fclose(LHS);

  FILE* RHS=fopen("RHSvector", "w");
  fprintf(RHS, "%d\n", N);
  for(int i=0; i<N; i++)
      fprintf(RHS, "%0.16e\n", b[i]);
  fclose(RHS);
  */

  
  double * rhs  = new double[N];
  double * res  = new double[N];
  int *permr = new int[N];
  int *rpermr = new int[N];
  int *permp = new int[2*N];

  for(i=0 ; i<N ; i++) 
    {
      permr[i] = rpermr[i] = permp[i+N] = i+1;
      res[i] = 0.0;
    }

  if(renumber == "rcmk")
    {
      a_rcmk    = new double[nbz];
      jptr_rcmk    = new int[N+1];
      ai_rcmk      = new int[nbz];
      int *mask         = new int[nbz];
      int *levels       = new int[N];
      i = j = k = 1;
      mat.cmkreord(N, a, ai, jptr, a_rcmk, ai_rcmk, jptr_rcmk, &i,
		   permr, mask, &j, &k, rpermr, levels);      
      double *vv = new double[nbz];
      mat.sort_col(N, a_rcmk, ai_rcmk, jptr_rcmk, mask, vv);
      delete [] vv;
      delete [] mask;
      delete [] levels;
      a = a_rcmk;
      jptr = jptr_rcmk;      
      ai = ai_rcmk;
    } 


  if(prec == "ilu")
    {
      int NbElm_ILU = nbz + N + 1;
      int* iw = (new int[ N + 1 ]);
      alu = new double [NbElm_ILU];
      jlu = new int [NbElm_ILU];
      ju = new  int [N+1];
      ilu0_ (&N,a,ai,jptr,alu,jlu,ju,iw,&ierr); 	       
      delete [] iw;
    }
  else if (prec == "ilut")
    {
      int nbfill = N;
      double droptol = 0.0;
      int NbElm_ILU = N*N*10;
   
      int* iw =  (new int   [ 2 * (N + 1) ]);
      double* dw = (new double[ (N + 1) ]);      
      alu = new double [NbElm_ILU];
      jlu = new int [NbElm_ILU];
      ju = new  int [N+1];

      ilut_ (&N,a,ai,jptr,&nbfill,&droptol,alu,jlu,ju,
	     &NbElm_ILU, dw, iw, &ierr);

      delete [] iw;
      delete [] dw;
    }
  /*
  else
    {
      throw 1;
    }
  */

  switch (ierr){
  case  0: break;
  case -1: fprintf(stderr, "# FATAL ERROR ON ILU : INPUT MATRIX MAY BE WRONG\n"); throw 1;
  case -2: fprintf(stderr, "# FATAL ERROR ON ILU : MATRIX L OVERFLOWS WORK ARRAY 'al'\n"); throw 1;
  case -3: fprintf(stderr, "# FATAL ERROR ON ILU : MATRIX U OVERFLOWS WORK ARRAY 'alu'\n"); throw 1;
  case -4: fprintf(stderr, "# FATAL ERROR ON ILU : ILLEGAL VALUE OF NB_FILL\n"); throw 1;
  case -5: fprintf(stderr, "# FATAL ERROR ON ILU : ZERO ROW ENCOUNTERED\n"); throw 1;
  default: fprintf(stderr, "# FATAL ERROR ON ILU : ZERO PIVOT ON LINE %d = %lf\n",ierr,
		   mat.GetMatrix(ierr,ierr));throw 1;
  }

  /* Pivoting rhs */    
    for(i=0;i<N;i++) {
      rhs[i] = b[rpermr[i] - 1];
    }
  
 
  
  /* Iterations */
  
  ipar[ 0] = 1;		// 1 for matvec with A, 2 with A^T, etc. see ipar(1) in iters.c 
  ipar[ 1] = 1;		// 1 => left preconditioner (2 => right, 0 => no preconditioner)
  ipar[ 4] = 100 ;	// Krylov size
  // changed
  ipar[ 5] = 10000; // nb iter max
  fpar[ 0] = 1.e-8; // stopping test
  fpar[ 1] = 1.e-8; // stopping test
  fpar[ 1] = 0.0;    
  fpar[10] = 0.0;

  if (accl == "gmres") {
	void sort_col(int n,double * a, int * ja, int *ia, int * iw, double * rw);

	int six = 0;
	ipar[3] = (((N + 1)+3) *(ipar[4]+2) + (ipar[4]+1)*ipar[4]/2);

//	if (prec == "noPC")	ipar[1] = 0;	// added by T. Bui 4/03

	//std::auto_ptr<double> work (new double[ ipar[3] ]);
	double* work = (new double[ ipar[3] ]);
	pgmres_ (&N, &ipar[4], rhs, res, work, &fpar[0],
			&ipar[5], &six, a, ai, jptr, 
			alu, jlu, ju, &ierr);
	if(ierr == 0) {
          int me;
          MPI_Comm_rank(MPI_COMM_WORLD,&me);
	  if(me == 0) {
	    std::cout << "iterative solver " << accl << " has converged\n";
  //	    std::cout << "iterative solver " << accl << " has converged to " << fpar[3] << " in " << ipar[6] << " iterations\n";
	  }
	}
	else if(ierr == 1) {
	  std::cout << "iterative solver " << accl << " has not converged in " << ipar[5] << " iterations\n"; 	  
	}
	else if(ierr == -1) {
	  std::cout << "iterative solver " << accl << " has faced a break-down\n"; 	  
	}      
        delete [] work;
  }
  else if (accl == "cg" || accl == "bcg") {
	int solver;
	if (accl == "cg") {
		ipar[3] = (N+1)*5;
		if (prec == "noPC")	ipar[1] = 0;	// added by T. Bui 4/03
		solver = 1;
	}
	else if (accl == "bcg") {
		ipar[3] = (N+1)*7;
		solver = 2;
	}

	//std::auto_ptr<double> work (new double[ ipar[3] ]);
	double* work = (new double[ ipar[3] ]);
	runrc_ (&N, rhs, res, ipar, fpar, work, res ,a, ai, jptr, alu, jlu, ju, &solver);
      
	if(ipar[0] == 0) {
	  std::cout << "iterative solver " << accl << " has converged to " << fpar[5] << " in " << ipar[6] << " iterations\n"; 	  
	}
	else if(ipar[0] == -1) {
	  std::cout << "iterative solver " << accl << " has not converged in " << ipar[6] << " iterations, residual is " << fpar[5] << "\n"; 	  
	}
	else if(ipar[0] == -3) {
	  std::cout << "iterative solver " << accl << " has faced a break-down\n"; 	  
	}      
        delete [] work;
  }  

// back solve
else if (accl == "backsolve") {
	lusol_ (&N, rhs, res, alu, jlu, ju);
  }
else {
  printf("ERROR: Unrecognized Sparskit Solver Option\n");
  throw 1;
}

  for (i=0;i<N;i++) {
	j = permr[i] - 1;
	k = permp[j+N] - 1;        
	x[i] = res[k];
  }
  
  delete [] rhs;
  delete [] res;
  delete [] alu;  
  delete [] jlu;
  delete [] ju;
  delete [] rpermr;
  delete [] permr;
  delete [] permp;    
  if (renumber == "rcmk") {
	delete [] a_rcmk;
	delete [] jptr_rcmk;
	delete [] ai_rcmk;
  }    
}
