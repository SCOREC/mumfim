#include "SparskitLinearSystemSolver.h"
#include "CSR_Vector.h"
#include "CSR_Matrix.h"


#include <assert.h>

namespace SCOREC {
  namespace Solver {  

    SparskitLinearSystemSolver::SparskitLinearSystemSolver(int nbRows) :
      theMatrix(NULL),
      theVector(NULL),
      theSolution(NULL)
    { }

    void SparskitLinearSystemSolver::LinearSystem_Initialize(int num_local_unkowns,
                                                             int num_global_unknowns,
							     int global_offset,
							     int * nonzero_per_row)
    {
      theMatrix = new CSR_Matrix(num_local_unkowns);
      theVector = new CSR_Vector(num_local_unkowns);
      theSolution = new CSR_Vector(num_local_unkowns);
    }

    void SparskitLinearSystemSolver::LinearSystem_AddToMatrix(int row, int col, double value)
    {
      // FORTRAN backend, so indexing starts at one
      theMatrix->AddMatrix(row+1,col+1,value);
    }
  
    void SparskitLinearSystemSolver::LinearSystem_AddToVector(int row, double value)
    {
      // same here
      theVector->AddVal(row+1,value);
    }
  
    bool SparskitLinearSystemSolver::LinearSystem_Zero( )
    {
      theMatrix->ZeroMatrix();
      theVector->ZeroArray();
      return true;
    }

    bool SparskitLinearSystemSolver::LinearSystem_ZeroMatrix( )
    {
      theMatrix->ZeroMatrix();
      return true;
    }
    bool SparskitLinearSystemSolver::LinearSystem_ZeroVector( )
    {
      theVector->ZeroArray();
      return true;
    }

    void SparskitLinearSystemSolver::LinearSystem_GetVector(double *& result)
    { 
      if (theVector)
	result = theVector->GetArray();
      else
	result = NULL;
    }

    void SparskitLinearSystemSolver::LinearSystem_SetVector(const double * x)
    { 
      //memcpy(&theVector[0],x,size*sizeof(double));
      if (theVector) {
	int size = theVector->size();
	for (int i = 0; i < size; i++)
	  (*theVector)[i] = x[i];
      }
    }

    void SparskitLinearSystemSolver::LinearSystem_Solve()
    {
      SPARSKIT_LINEAR_SOLVER_("rcmk",
			      "ilut",
			      "gmres",
			      *theMatrix, *theVector, *theSolution);
    }
  
    void SparskitLinearSystemSolver::LinearSystem_GetSolution(double *& sol)
    {
      if(theMatrix)
	sol = theSolution->GetArray();
      else
	sol = NULL;
    }

    void SparskitLinearSystemSolver::LinearSystem_PrintMatrix(std::ostream & outstream) 
    {
      if(theMatrix)
	theMatrix->OutputMatrixMatlabFormat(outstream);  
    }
  
    void SparskitLinearSystemSolver::LinearSystem_PrintVector(std::ostream & outstream) 
    {
      if(theVector)
	theVector->PrintArray(outstream);
    }

    SparskitLinearSystemSolver::~SparskitLinearSystemSolver( )
    {
      if (theMatrix)
      {
	delete theMatrix;
	delete theVector;
	delete theSolution;
      }
    }
  
  } 
}//end of namespace SCOREC_Solver
