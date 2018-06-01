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
   * @note Sparskit is a single-precision solver, so don't use this for anything where you're
   *  trying to converge past 10e-8
   */
  LasSolve * createSparskitLUSolve(SparskitBuffers * b, double eps = 0.0);
  LasSolve * createSparskitQuickLUSolve(SparskitBuffers * b);
  LasSolve * createSparskitQuickLUSolve(LasSolve * slvr);
  /*
   * sparse defines wheter to print the matrix in sparse mode, or as full matrix
   */
  void printSparskitMat(std::ostream &, Mat * m, bool sparse=false);
  void printSparskitVec(std::ostream &, Vec * v);
  double getSparskitMatValue(Mat *, int rr, int cc);
  void setSparskitMatValue(Mat *, int rr, int cc, double vl);
  double * getskMatNNZArray(Mat * m);
  SparskitBuffers * getSparskitBuffers(LasSolve * slv);
}
#endif
