#ifndef BIO_CSR_H_
#define BIO_CSR_H_

/**
 * @todo Bill : move to amsi
 */

#include <apfNumbering.h>
#include <vector>

namespace bio
{
  /**
   * A utility class to describe the sparse structure of a
   *  compressed sparse row (CSR) format matrix,
   *  can be used on in conjunction with any linear storage container
   *  to index a sparse matrix.
   */
  class CSR
  {
  private:
    int neq;
    int nnz;
    std::vector<int> rws;
    std::vector<int> cls;
    CSR();
  public:
    CSR(int ne, int nnz, int * rs, int * cs);
    int getNumEqs() const { return neq; }
    int getNumNonzero() const { return nnz; }
    int operator()(int,int) const;
    int * getRows() { return &rws[0]; }
    int * getCols() { return &cls[0]; }
  };

  /**
   * Construct a CSR sparse matrix structure to manage a sparse matrix based
   *  on the dof numbering of an apf field.
   */
  CSR * createCSR(apf::Numbering * num, int ndofs);

  /**
   *  Produce a full version of the sparse matrix, this is typically
   *   only used for debugging purposes.
   * @param csr A csr object describing the structure of the sprs_mat
   * @param sprs_mat A buffer of length csr->getNumNonzeros() containing the matrix values
   * @param fll_mat A preallocated buffer of length csr->getNumEqs()^2 which will contain
   *                all nonzero and zero values of the matrix;
   */
  void constructFullMatrix(CSR * csr,double * sprs_mat,double * fll_mat);
}
#endif
