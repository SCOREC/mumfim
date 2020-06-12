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
    constexpr double h = 1E-5;
    double D1 = sqrt(1.0 / (1 - 2 * h));
    constexpr double l1 = 1;
    constexpr double l2 = (1 + h) / (1 - h * h);
    constexpr double l3 = (1 - h) / (1 - h * h);
    double l2pl3 = 0.5 * (sqrt(l2) + sqrt(l3));
    double l2ml3 = 0.5 * (sqrt(l2) - sqrt(l3));
    // assume 3D
    for(int i=0; i<6; ++i)
    {
      // compute V from V = sqrt((I-2e)^-1)
      // symmetric e=1/2(e+e^T)
      DeformationGradient Fappd;
      switch (i)
      {
        case 0:
          Fappd = DeformationGradient(D1, 0, 0, 0, 1, 0, 0, 0, 1);
          break;
        case 1:
          Fappd = DeformationGradient(1, 0, 0, 0, D1, 0, 0, 0, 1);
          break;
        case 2:
          Fappd = DeformationGradient(1, 0, 0, 0, 1, 0, 0, 0, D1);
          break;
        case 3:
          Fappd =
              DeformationGradient(l1, 0, 0, 0, l2pl3, l2ml3, 0, l2ml3, l2pl3);
          break;
        case 4:
          Fappd =
              DeformationGradient(l2pl3, 0, l2ml3, 0, l1, 0, l2ml3, 0, l2pl3);
          break;
        case 5:
          Fappd =
              DeformationGradient(l2pl3, l2ml3, 0, l2ml3, l2pl3, 0, 0, 0, l1);
          break;
      }
      run(Fappd, sigma1, false);
      for(int j=0; j<6; ++j)
      {
        // j is the row, i is the column of 6x6 C
        C[j*6+i] = (sigma1[j]-curStress[j])/h;
      }
    }
  }
}
