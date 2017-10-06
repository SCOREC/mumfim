#include "lasSparskit.h"
#include "lasSparskitExterns.h"
#include "lasCSR.h"
#include <cassert>
#include <cstring>
#include <iostream>
namespace las
{
  class SparskitLU : public LasSolve
  {
  private:
    SparskitBuffers * bfrs;
  public:
    SparskitLU(SparskitBuffers * b) : bfrs(b) {}
    virtual void solve(Mat * k, Vec * u, Vec * f);
  };
  class SparskitOps : public LasOps
  {
  public:
    virtual void zero(Mat * m);
    virtual void zero(Vec * v);
    virtual void assemble(Vec * v, int cnt, int * rws, double * vls);
    virtual void assemble(Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls);
    virtual void set(Vec * v, int cnt, int * rws, double * vls);
    virtual void set(Mat * m, int cntr, int * rws, int cnts, int * cls, double * vls);
    virtual double norm(Vec * v);
    virtual double dot(Vec * v0, Vec * v1);
    virtual void axpy(double a, Vec * x, Vec * y);
    virtual void get(Vec * v, double *& vls);
    virtual void restore(Vec * , double *& ) {};
  };
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
    ~skVec()
    {
      delete [] vls;
    }
    double & operator[](int idx)
    {
      assert(idx >= 0 && idx < cnt);
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
    double * vls;
    CSR * csr;
    double zero;
  public:
    skMat(CSR * c)
      : vls(new double [c->getNumNonzero()]())
      , csr(c)
      , zero(0.0)
    { }
    ~skMat()
    {
      delete [] vls;
    }
    double & operator()(int rr, int cc)
    {
      zero = 0.0;
      if(rr < 0 || cc < 0)
        return zero;
      int idx = (*csr)(rr,cc);
      if(idx < 0)
        return zero;
      return vls[idx];
    }
    CSR * getCSR()
    {
      return csr;
    }
  };
  skMat * getSparskitMatrix(Mat * m)
  {
    return reinterpret_cast<skMat*>(m);
  }
  skVec * getSparskitVector(Vec * v)
  {
    return reinterpret_cast<skVec*>(v);
  }
  Mat * createSparskitMatrix(CSR * csr)
  {
    return reinterpret_cast<Mat*>(new skMat(csr));
  }
  void deleteSparskitMatrix(Mat * m)
  {
    skMat * mat = getSparskitMatrix(m);
    delete mat;
  }
  Vec * createSparskitVector(int n)
  {
    return reinterpret_cast<Vec*>(new skVec(n));
  }
  void deleteSparskitVector(Vec * v)
  {
    skVec * vec = getSparskitVector(v);
    delete vec;
  }
  LasOps * initSparskitOps()
  {
    static SparskitOps * ops = NULL;
    if(ops == NULL)
      ops = new SparskitOps;
    return ops;
  }
  LasSolve * createSparskitLUSolve(SparskitBuffers * b)
  {
    return new SparskitLU(b);
  }
  void printSparskitMat(std::ostream & o, Mat * mi)
  {
    skMat * m = getSparskitMatrix(mi);
    int ndofs = m->getCSR()->getNumEqs();
    for(int rr = 0; rr < ndofs; ++rr)
    {
      for(int cc = 0; cc < ndofs; ++cc)
      {
        o << (*m)(rr,cc) << ' ';
      }
      o << '\b' << std::endl;
    }
  }
  // CLASS MEMBER FUNCTION DEFINITIONS
  void SparskitLU::solve(Mat * k, Vec * u, Vec * f)
  {
    bfrs->zero();
    double tol = 1e-6;
    skMat * mat = getSparskitMatrix(k);
    skVec * uv = getSparskitVector(u);
    skVec * fv = getSparskitVector(f);
    CSR * csr = mat->getCSR();
    int bfr_lng = bfrs->matrixLength();
    int ierr = 0;
    int ndofs = csr->getNumEqs();
    ilut_(&bfr_lng,
          &(*mat)(0,0),
          csr->getCols(),
          csr->getRows(),
          &ndofs,
          &tol,
          bfrs->matrixBuffer(),
          bfrs->colsBuffer(),
          bfrs->rowsBuffer(),
          &bfr_lng,
          bfrs->doubleWorkBuffer(),
          bfrs->intWorkBuffer(),
          &ierr);
    if(ierr != 0)
    {
      std::cerr << "ERROR: ilut_ returned error code " << ierr << std::endl;
      return;
    }
    lusol_(&ndofs,
           &(*fv)[0],
           &(*uv)[0],
           bfrs->matrixBuffer(),
           bfrs->colsBuffer(),
           bfrs->rowsBuffer());
  }
  void SparskitOps::zero(Mat * m)
  {
    skMat * mat = getSparskitMatrix(m);
    memset(&(*mat)(0,0),0.0,mat->getCSR()->getNumNonzero()*sizeof(double));
  }
  void SparskitOps::zero(Vec * v)
  {
    skVec * vec = getSparskitVector(v);
    memset(&vec[0],0.0,vec->size()*sizeof(double));
  }
  void SparskitOps::assemble(Vec * v, int cnt, int * rws, double * vls)
  {
    skVec * vec = getSparskitVector(v);
    for(int ii = 0; ii < cnt; ++ii)
      (*vec)[rws[ii]] += vls[ii];
  }
  void SparskitOps::assemble(Mat * m, int cntr, int * rws, int cntc, int * cols, double * vls)
  {
    skMat * mat = getSparskitMatrix(m);
    for(int ii = 0; ii < cntr; ++ii)
      for(int jj = 0; jj < cntc; ++jj)
      {
        double vl = vls[ii * cntc + jj];
        if(vl != 0.0) // don't want to attempt to access zero locations in a CSR matrix
          (*mat)(rws[ii],cols[jj]) += vls[ii * cntc + jj];
      }
  }
  void SparskitOps::set(Vec * v, int cnt, int * rws, double * vls)
  {
    skVec * vec = getSparskitVector(v);
    for(int ii = 0; ii < cnt; ++ii)
      (*vec)[rws[ii]] = vls[ii];
  }
  void SparskitOps::set(Mat * m, int cntr, int * rws, int cntc, int * cols, double * vls)
  {
    skMat * mat = getSparskitMatrix(m);
    for(int ii = 0; ii < cntr; ++ii)
      for(int jj = 0; jj < cntc; ++jj)
        (*mat)(rws[ii],cols[jj]) = vls[ii * cntc + jj];
  }
  double SparskitOps::norm(Vec * v)
  {
    skVec * vec = getSparskitVector(v);
    double nrm = 0.0;
    for(int ii = 0; ii < vec->size(); ++ii)
      nrm += (*vec)[ii] * (*vec)[ii];
    nrm = sqrt(nrm);
    return nrm;
  }
  double SparskitOps::dot(Vec * v0, Vec * v1)
  {
    skVec * vec0 = getSparskitVector(v0);
    skVec * vec1 = getSparskitVector(v1);
    int sz0 = vec0->size();
    assert(sz0 == vec1->size());
    double dt = 0.0;
    for(int ii = 0; ii < sz0; ++ii)
      dt += (*vec0)[ii] * (*vec1)[ii];
    return dt;
  }
  // y = ax + y
  void SparskitOps::axpy(double a, Vec * x, Vec * y)
  {
    skVec * vx = getSparskitVector(x);
    skVec * vy = getSparskitVector(y);
    int szx = vx->size();
    assert(szx == vy->size());
    for(int ii = 0; ii < szx; ++ii)
      (*vy)[ii] = a * (*vx)[ii] + (*vy)[ii];
  }
  void SparskitOps::get(Vec * v, double *& vls)
  {
    skVec * vec = getSparskitVector(v);
    vls = &(*vec)[0];
  }
}
