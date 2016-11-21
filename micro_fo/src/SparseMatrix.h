#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_
#include "FiberNetwork.h"
//#include <memory>
#include <vector>
namespace bio
{
  class SparseMatrix
  {
  private:
    int num_nonzeros;
    bool finalized;
    std::vector<int> rows;
    std::vector<int> cols;
    SparseMatrix();
    int findSpot(int,int) const;
  public:
    SparseMatrix(int);
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
    void debugFullMatrix(double * matrix, std::vector<double> & fmatrix, int n);
    // below here will break until finalize is called
    //TODO: protect against that..
    int sparseLocation(int,int) const;
    int * getRows() {return &rows[0];}
    int * getCols() {return &cols[0];}
  };
}
#endif
