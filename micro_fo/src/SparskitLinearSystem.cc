#include "SparskitLinearSystem.h"
#include "CSR.h"

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

namespace bio
{
  skVec makeVec(int rws) { return new double[rws]; }
  void destroyVec(skVec vec) { delete [] vec; }
  
  skMat::skMat(CSR * ss)
    : mat(NULL)
    , sttr(ss)
  {
    assert(ss);
    mat = new double[sttr->getNumNonzero()];
  }

  skMat::~skMat()
  {
    delete [] mat;
  }
  double skMat::operator()(int r, int c) const
  {
    return mat[(*sttr)(r,c)];
  }
  double & skMat::operator()(int r,int c)
  {
    return mat[(*sttr)(r,c)];
  }

  void skSolver::solve(skMat k, skVec u, skVec f)
  {
    CSR * ss = k.getSparseStructure();
    /*
    ilut_(&ndofs,
	  &k,
	  ss->getCols(),
	  ss->getRows(),
	  &ndofs,
	  &tol,
	  matrixbuffer,
	  colbuffer,
	  rowbuffer,
	  buffer_length,
	  doubleworkbuffer,
	  intworkbuffer,
	  &ierr);
    if(ierr)
      std::cerr << "ilut_ error code: " << ierr << std::endl;
    
    lusol(&ndofs,
	  f,
	  u,
	  matrixbuffer,
	  colsbuffer,
	  rowsbuffer);
    */
  }
}
