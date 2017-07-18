#ifndef LAS_SPARSKIT_H_
#define LAS_SPARKSIT_H_
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
  LasSolve * createSparskitLUSolve(SparskitBuffers * b);
}
#endif
