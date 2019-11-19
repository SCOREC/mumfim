#include "bioNeoHookeanRVEAnalysis.h"
#include <apfMatrixUtil.h>
namespace bio
{
  bool NeoHookeanRVEAnalysis::run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords)
  {
    if (!update_coords)
      F_old = F;
    F = apf::Matrix3x3(dfmGrd[0], dfmGrd[1], dfmGrd[2],
                  dfmGrd[3], dfmGrd[4], dfmGrd[5],
                  dfmGrd[6], dfmGrd[7], dfmGrd[8]);
    // compute the left cauchy green deformation tensor
    leftCauchyGreen = apf::DynamicMatrix(3, 3);   // leftCauchyGreen.zero();
    apf::DynamicMatrix FT(3, 3);
    FT.zero();
    apf::transpose(fromMatrix(F), FT);
    apf::multiply(fromMatrix(F), FT, leftCauchyGreen);
    computeCauchyStress(sigma);
    if(!update_coords)
      F = F_old;
    else
    {
      for(int i=0; i<6; ++i)
        curStress[i] = sigma[i];
    }
    return true;
  }
  void NeoHookeanRVEAnalysis::computeCauchyStress(double sigma[6])
  {
    detF = apf::getDeterminant(F);
    apf::Matrix3x3 Cauchy;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Cauchy[i][j] = ShearModulus / detF * (leftCauchyGreen(i, j) - (i == j)) +
                       lambda / detF * log(detF) * (i == j);
      }
    }
    amsi::mat2VoigtVec(3, Cauchy, sigma);
  }
  void NeoHookeanRVEAnalysis::computeMaterialStiffness(double C[36])
  {
    detF = apf::getDeterminant(F);
    double lambda_prime = lambda / detF;
    double mu_prime = (ShearModulus - lambda * log(detF)) / detF;
    apf::DynamicMatrix D(6, 6);
    D.zero();
    D(0, 0) = lambda_prime + (2.0 * mu_prime);
    D(0, 1) = lambda_prime;
    D(0, 2) = lambda_prime;
    D(1, 0) = lambda_prime;
    D(1, 1) = lambda_prime + (2.0 * mu_prime);
    D(1, 2) = lambda_prime;
    D(2, 0) = lambda_prime;
    D(2, 1) = lambda_prime;
    D(2, 2) = lambda_prime + (2.0 * mu_prime);
    D(3, 3) = mu_prime;
    D(4, 4) = mu_prime;
    D(5, 5) = mu_prime;
    double Ctest[36];
    for(int i=0; i<6; ++i)
      for(int j=0; j<6; ++j)
        C[i*6+j]= D(i,j);
        //C[i*6+j]= D(i,j);
   //RVEAnalysis::computeMaterialStiffness(C);
   //for(int i=0; i<36; ++i)
   //  std::cout<<fabs(Ctest[i] -C[i])<<" ";
   //std::cout<<std::endl;
   
  }
  NeoHookeanRVEAnalysis * initNeoHookeanRVEAnalysisFromMultiscale(micro_fo_params & prm)
  {
    double youngs_modulus = prm.data[YOUNGS_MODULUS];
    double poisson_ratio = prm.data[NONLINEAR_PARAM];
    if(poisson_ratio > 0.5)
    {
      std::cerr<<"Poisson ratio must be less than 0.5"<<std::endl;
      MPI_Abort(AMSI_COMM_SCALE, 1);
    }
    return new NeoHookeanRVEAnalysis(youngs_modulus, poisson_ratio);
  }
}
