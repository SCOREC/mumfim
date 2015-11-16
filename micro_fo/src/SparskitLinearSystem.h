#ifndef BIO_SPARSKIT_LINEAR_SYSTEM_H_
#define BIO_SPARSKIT_LINEAR_SYSTEM_H_

#include "Sparskit_Externs.h"

namespace bio
{

  Vector * castVector(double *);
  Matrix * castMatrix(double *);

  class SparskitOps : public LSOps
  {
    
  };

  class SparskitLinearSystem : public LinearSystem
  {
  private:

  public:
    Vector * getVector();
    Vector * getSolution();
    Matrix * getMatrix();
    LSOps* getOps();
  };

  class SparskitSolver : public LinearSolver
  {
  };

  
}

#endif
