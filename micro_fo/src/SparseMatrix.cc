#include "SparseMatrix.h"
#include <cstdlib>
#include <iostream>
namespace bio
{
  SparseMatrix::SparseMatrix(int num_eqs)
    : num_nonzeros(0)
    , finalized(false)
    , rows(num_eqs+1)
    , cols(num_eqs*10)
  {
    rows.assign(num_eqs+1,1);
    cols.assign(num_eqs*sqrt(num_eqs),0);
  }
  bool SparseMatrix::needMtxElem(int col,int row)
  {
    bool result = false;
    if(!finalized)
    {
      result = true;
      for(int ii = rows.at(row); ii < rows.at(row+1); ii++)
      {
        if(cols.at(ii-1) == (col+1))
        {
          result = false;
          break;
        }
      }
    }
    return result;
  }
  void SparseMatrix::addMtxElem(int col, int row,int n)
  {
    if(!finalized)
    {
      int ii = findSpot(col,row) - 1;
      for(int jj = row + 1; jj <= n; jj++)
        rows.at(jj)++;
      for(int jj = rows.at(n) - 2; jj > ii; jj--)
        cols.at(jj) = cols.at(jj - 1);
      cols.at(ii) = col + 1;
      num_nonzeros++;
    }
  }
  int SparseMatrix::findSpot(int col, int row) const
  {
    int ii;
    for(ii = rows.at(row); (ii < rows.at(row+1)) && (cols.at(ii-1) < (col + 1)); ii++) {}
    return ii;
  }
  void SparseMatrix::finalize()
  {
    cols.resize(num_nonzeros); // only store exactly as much as is needed
    finalized = true;
  }
  int SparseMatrix::sparseLocation(int col, int row) const
  {
    int result = -1;
    int kk;
    for(kk = rows.at(row); (kk < rows.at(row+1)) && (cols.at(kk-1) < (col+1)); kk++){}
    if(cols.at(kk - 1) == (col + 1))
      result = kk - 1;
    else
    {
      std::cerr << "ERROR: Lost data point " << row << "," << col << std::endl;
      std::exit(EXIT_FAILURE);
    }
    return result;
  }
  // Just for debugging - put sparse matrix representation into vector in full representation
  void SparseMatrix::debugFullMatrix(double * matrix, std::vector<double> & fmatrix, int n)
  {
    if(finalized)
    {
      finalized = false; // Just to be able to use needMtxElem
      fmatrix.resize(n*n);
      for(int ii=0;ii<n;ii++)
        for(int jj=0;jj<n;jj++)
        {
          if(needMtxElem(jj,ii))
            fmatrix[ii*n + jj] = 0;
          else
            fmatrix[ii*n + jj] = matrix[sparseLocation(jj,ii)];
        }
      finalized = true;
    }
  }
}
