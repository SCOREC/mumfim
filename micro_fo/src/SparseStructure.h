#ifndef BIO_SPARSE_STRUCTURE_H_
#define BIO_SPARSE_STRUCTURE_H_

#include "FiberNetwork.h"
#include <vector>

namespace bio
{
  class SparseStructure {
  private:
    int num_nonzeros;
    bool finalized;
    std::vector<int> rows;
    std::vector<int> cols;
    int findSpot(int,int) const;
    SparseStructure();
  public:
    SparseStructure(int);
    int numNonzeros() const
    {
      if(finalized)
	return num_nonzeros;
      else
	return -1;
    }
    bool needMtxElem(int,int);
    void addMtxElem(int,int,int);
    void finalize();
    bool isFinal() {return finalized;}
    int sparseLocation(int,int) const;
    int * getRows()
    {
      int * result = NULL;
      if(finalized)
	result = &rows[0];
      return result;
    }
    int * getCols()
    {
      int * result = NULL;
      if(finalized)
	result = &cols[0];
      return result;
    }
  };

  void extractFullMatrix(SparseStructure & ss,
			 double * sprs_mtrx,
			 double * fll_mtrx,
			 int n);

}

#endif
