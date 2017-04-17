#include "lasSparskit.h"
#include "lasCSR.h"
namespace las
{
  class skVec
  {
  private:
    double * vls;
    int cnt;
  public:
    skVec(int n)
      : vls(new double [n])
      , cnt(n)
      { }
    double & operator[](int idx)
    {
      assert(idx >= 0 && idx < n);
      return vls[idx];
    }
    int size()
    {
      return cnt;
    }
  };
  class skMat
  {
  private:
  };
  Mat * createSparskitMatrix(CSR * csr)
  {
    return reinterpret_cast<Mat*>(new skMat(csr));
  }
  Vec * createSparskitVector(int n)
  {
    return reinterpret_cast<Vec*>(new skVec(n));
  }
  LasOps * initSparskitOps()
  {
    static SparskitOps * ops = NULL;
    if(ops == NULL)
      ops = new SparskitOps;
    return ops;
  }
  CSRMat * getSparskitMatrix(Mat * m)
  {
    return reinterpret_cast<skMat*>(m);
  }
  CSRMat * getSparskitVector(Vec * v)
  {
    return reinterpret_cast<skVec*>(v);
  }
  void SparskitOps::zero(Mat * m)
  {

  }
  void SparskitOps::zero(Vec * v)
  {
    skVec * vec = getSparskitVector()
    memset(&v[0],0.0,v->size());
  }
  void SparskitOps::assemble(Vec * v, int cnt, int * rws, double * vls)
  {
    for(int ii = 0; ii < cnt; ++ii)
      v[rws[ii]] += vls[ii];
  }
  void SparskitOps::assemble(Mat * m, int cntr, int * rws, int cntc, int * cols, double * vls)
  {

  }
  void SparskitOps::set(Vec * v, int cnt, int * rws, double * vls)
  {
    for(int ii = 0; ii < cnt; ++ii)
      v[rws[ii]] = vls[ii];
  }
  void SparskitOps::set(Mat * m, int cnts, int * rws, int cntc, int * cols, double * vls)
  {

  }
  void SparskitOps::solve(Mat * k, Vec * u, Vec * f)
  {

  }
  double SparskitOps::norm(Vec * v)
  {

  }
  double SparskitOps::dot(Vec * v0, Vec * v1)
  {

  }
  // y = ax + y
  void SparskitOps::axpy(double a, Vec * x, Vec * y)
  {

  }
  void SparskitOps::get(Vec * v, double *& vls)
  {
    vls = v;
  }
}
