#include "lasCSR.h"
#include <cstring>
#include <iostream>
namespace las
{
  CSR::CSR(int ne, int nz, int * rs, int * cs)
    : neq(ne)
    , nnz(nz)
    , rws(rs,rs+ne+1)
    , cls(cs,cs+nz)
  { }
  int CSR::operator()(int rw, int cl) const
  {
    int result = -1;
    int fst = rws[rw] - 1; // 1-indexing    v also here
    while((fst < rws[rw+1] - 2) && (cls[fst] - 1 < cl))
      ++fst;
    if(cls[fst] - 1 == cl) // 1-indexing
      result = fst;
    else
      result = -1;
    return result;
  }
  void constructFullMatrix(CSR * csr,double * sprs_mat,double * fll_mat)
  {
    int neq = csr->getNumEqs();
    memset(&fll_mat[0],0,neq*neq*sizeof(double));
    int lc = -1;
    for(int ii = 0; ii < neq; ii++)
      for(int jj = 0; jj < neq; jj++)
        if((lc = (*csr)(ii,jj)) != -1)
          fll_mat[ii*neq + jj] = sprs_mat[lc];
  }
}
