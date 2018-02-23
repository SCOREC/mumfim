#ifndef LAS_SPARSKIT_H_
#define LAS_SPARSKIT_H_
#include "lasCSR.h"
#include "lasSparskitExterns.h"
#include <amsiLAS2.h>
namespace las
{
  Mat * createSparskitMatrix(CSR * csr);
  void deleteSparskitMatrix(Mat * m);
  Vec * createSparskitVector(int n);
  void deleteSparskitVector(Vec * v);
  LasOps * initSparskitOps();
  /**
   * @note Sparskit is a single-precision solver, so don't use this for anything nonlinear where you're
   *  trying to converge past 10e-8
   */
  LasSolve * createSparskitLUSolve(SparskitBuffers * b);
  LasSolve * createSparskitQuickLUSolve(SparskitBuffers * b);
  LasSolve * createSparskitQuickLUSolve(LasSolve * slvr);
  void printSparskitMat(std::ostream &, Mat * m);
  double getSparskitMatValue(Mat *, int rr, int cc);
  void setSparskitMatValue(Mat *, int rr, int cc, double vl);
}
#endif
