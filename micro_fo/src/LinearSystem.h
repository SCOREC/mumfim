#ifndef BIO_LINEAR_SYSTEM_H_
#define BIO_LINEAR_SYSTEM_H_

#include <apfDynamicMatrix.h>
#include <apfNew.h>

namespace bio
{
  class Vector;
  class Matrix;

  //Vector * castVector(); specialize for different linear solver backends...
  //Matrix * castMatrix();

  /**
   * Performs general operations on the objects of a linear system, vectors and matrices
   */
  class LSOps
  {
  public:
    virtual void setVector(Vector * vec,
			   const apf::DynamicVector & fe,
			   const apf::NewArray<int> & dofs,
			   int nedofs) = 0;
    virtual void addToVector(Vector * vec,
			     const apf::DynamicVector & fe,
			     const apf::NewArray<int> & dofs,
			     int nedofs) = 0;
    virtual void addToMatrix(Matrix * mat,
			     const apf::DynamicMatrix & ke,
			     const apf::NewArray<int> & dofs,
			     int nedofs) = 0;
    virtual double norm(Vector * vec) = 0;
  };

  class LinearSystem;
  /**
   * Handles solving a linear sytem
   */
  class LinearSolver
  {
  public:
    virtual void solve(LinearSystem *) = 0;
  };

  /**
   * Manages the data structures associated with a single linear system \f$ A x = b \f$.
   */
  class LinearSystem
  {
  protected:
  public:
    virtual Vector * getVector() = 0;
    virtual Vector * getSolution() = 0;
    virtual Matrix * getMatrix() = 0;
    virtual LSOps * getOps() = 0;
    virtual LinearSolver * getSolver() = 0;
  };
}

#endif
