/****************************************************************************** 

  (c) 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#ifndef H_SparskitLinearSystemSolver
#define H_SparskitLinearSystemSolver

#include "LAS.h"

#include <ostream>

class CSR_Matrix;
class CSR_Vector;

namespace SCOREC {
  namespace Solver {
  
    class SparskitLinearSystemSolver : public LinearSystemSolver 
    {
      
      CSR_Matrix * theMatrix;
      CSR_Vector * theVector;
      CSR_Vector * theSolution;
      
    public:
      virtual void LinearSystem_Initialize(int, int, int, int*);

      /// add an element val at given row and column in the system
      virtual void LinearSystem_AddToMatrix(int row, int col, double value );

      /// add an element in the right hand side
      virtual void LinearSystem_AddToVector(int row, double value );

      /// solve the system
      virtual void LinearSystem_Solve();
      
      /// zero the system
      virtual bool LinearSystem_Zero(); 
      
      /// retrieve the solution
      virtual void LinearSystem_GetSolution(double *&);

      virtual bool LinearSystem_ZeroMatrix();
      virtual bool LinearSystem_ZeroVector();
      
      virtual void LinearSystem_GetVector(double *&);
      virtual void LinearSystem_SetVector(const double *);

      virtual void LinearSystem_GetVectorNorm(double &);

      SparskitLinearSystemSolver ( int nbRows );
      ~SparskitLinearSystemSolver ( );

      virtual bool LinearSystem_SeparateSolve() {return false;}

      virtual void LinearSystem_PrintMatrix(std::ostream &);
      virtual void LinearSystem_PrintVector(std::ostream &);
    
  private:

      // no copy/assignment 
      SparskitLinearSystemSolver ( const SparskitLinearSystemSolver & )
      {};
      SparskitLinearSystemSolver operator = ( const SparskitLinearSystemSolver &)
	{return *this;}        
    };

  } 
} // end of namespace SCOREC::Solver

#endif
