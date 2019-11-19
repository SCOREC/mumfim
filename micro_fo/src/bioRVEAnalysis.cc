#include "bioRVEAnalysis.h"
#include "bioMicroFOParams.h"
#include <apfMatrix.h>
namespace bio
{
  RVEAnalysis::RVEAnalysis()
  {
    std::fill_n(&curStress[0], 6, 0);
  } 
  RVEAnalysis::RVEAnalysis(const RVEAnalysis & an)
  {
    for(int i=0; i<6; ++i)
      curStress[i] = an.curStress[i];
  }
  void RVEAnalysis::computeMaterialStiffness(double C[36])
  {
    std::vector<std::vector<int>> idx = {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};
    double sigma1[6];
    //double h = 1E-8;
    double h = 1E-5;
    // assume 3D
    for(int i=0; i<6; ++i)
    {
      // compute V from V = sqrt((I-2e)^-1)
      // compute I-2e
      apf::Matrix3x3 im2e(1,0,0,0,1,0,0,0,1);
      // symmetric e=1/2(e+e^T)
      im2e[idx[i][0]][idx[i][1]] = im2e[idx[i][1]][idx[i][0]] = (idx[i][0]==idx[i][1]) ? 1-2*h : -h;
      apf::Matrix3x3 im2einv = apf::invert(im2e);
      apf::Matrix3x3 V;
      apf::applyMatrixFunc(im2einv, &sqrt, V);
      DeformationGradient Fappd(V[0][0], V[0][1], V[0][2],
                                V[1][0], V[1][1], V[1][2],
                                V[2][0], V[2][1], V[2][2]);
      run(Fappd, sigma1, false);
      for(int j=0; j<6; ++j)
      {
        // j is the row, i is the column of 6x6 C
        C[j*6+i] = (sigma1[j]-curStress[j])/h;
      }
    }
  }
}
