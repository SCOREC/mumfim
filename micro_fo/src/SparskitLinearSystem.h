#ifndef BIO_SPARSKIT_LINEAR_SYSTEM_H_
#define BIO_SPARSKIT_LINEAR_SYSTEM_H_

#include <cassert>

namespace bio
{
  class CSR;
  
  typedef double* skVec;
  skVec makeVec(int rws);
  void destroyVec(skVec vec);

  class skMat
  {
  private:
    double * mat;
    CSR * sttr;
  public:
    skMat(CSR * ss);
    ~skMat();
    double operator()(int r, int c) const;
    double & operator()(int r,int c);
    CSR * getSparseStructure() const { return sttr; }
  };

  template <typename T1, typename T2>
    void setVecValues(skVec * vec, T1 & fe, T2 & rws, int nmrws, bool add = false);
  template <typename T1, typename T2>
    void setMatValues(skMat * mat, T1 & ke, T2 & rws, int nmrws, T2 & cls, int nmcls, bool add = false);

  class skSolver
  {
  public:
    void solve(skMat k, skVec u, skVec f);
  };
}
#include "SparskitLinearSystem_impl.h"
#endif

