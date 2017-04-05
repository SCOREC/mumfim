#ifndef LAS_SPARSKIT_H_
#define LAS_SPARKSIT_H_
#include <amsiLAS2.h>
namespace las
{
  Mat * createSparskitMatrix(int n);
  Vec * createSparskitVector(int n);
  LasOps * initSparskitOps();
  class SparskitOps
  {
  public:
    virtual void zero(Mat * m) = 0;
    virtual void zero(Vec * v) = 0;
    virtual void add(Vec * v, int cnt, int * rws, double * vls) = 0;
    virtual void add(Mat * m, int cntr, int * rws, int cntc, int * cols, double * vls) = 0;
    virtual void solve(const Mat * k, Vec * u, const Vec * f) = 0;
    virtual void dot(Vec * v0, Vec * v1) = 0;
    virtual void norm(Vec * v) = 0;
  };
}
#endif
