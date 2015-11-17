#ifndef BIO_SPARSE_STRUCTURE_H_
#define BIO_SPARSE_STRUCTURE_H_

#include <apfField.h>
#include <vector>

namespace bio
{

  SparseStructure * makeStructure(apf::Numbering * nm, int ndofs)
  {
    BuildCSRStruct mkr(nm,ndofs);
    mkr.apply(apf::getField(nm));
    SparseStructure * ss = mkr.finalize();
    return ss;
  }

  void extractFullMatrix(SparseStructure & ss,
			 double * sprs_mtrx,
			 double * fll_mtrx,
			 int n);

}

#endif
