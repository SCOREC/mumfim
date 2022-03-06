#ifndef MUMFIM_NEOHOOKEAN_INTEGRATOR_H_
#define MUMFIM_NEOHOOKEAN_INTEGRATOR_H_
#include <ElementalSystem.h>
#include <apfShape.h>
#include <math.h>  //natural log
#include <cstring>
#include "bioNonlinearTissue.h"
namespace mumfim
{
  // Parameters from V. Lai et al. Journal of Biomechanical Engineering, Vol
  // 135, 071007 (2013)
  //
  // The Voigt notation in this integrator are NOT consistent
  // (in different portions of the integrator)!
  // Fortunately, we are OK because this is an isotropic material,
  // but for consistency (with the rest of the code) 
  // we should use 1->11, 2->22, 3->33, 4->23, 5->13, 6->12
  class NeoHookeanIntegrator : public amsi::ElementalSystem {
    public:
    NeoHookeanIntegrator(NonlinearTissue *n,
			 apf::Field *displacements,
                         apf::Field *dfm_grd,
			 apf::Field *current_coords,
                         double youngs_modulus,
			 double poisson_ratio,
			 int o)
        : ElementalSystem(displacements, o)
        , current_integration_point(0)
        , analysis(n)
        , dim(0)
        , dfm_grd_fld(dfm_grd)
        , current_coords(current_coords)
        , ShearModulus(0.0)
        , PoissonsRatio(poisson_ratio)
    {
      ShearModulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    }
    // ce is the mesh element that corresponds to the current coordinates
    void inElement(apf::MeshElement *me)
    {
      ElementalSystem::inElement(me);
      msh = apf::getMesh(f);
      dim = apf::getDimension(me);
      current_integration_point = 0;
      // testing of new apf work
      // create a mesh element that tracks with current coords
      ccme = apf::createMeshElement(current_coords, apf::getMeshEntity(me));
      // current coordinate field w.r.t. current coordinate mesh elements
      // elements
      cccce = apf::createElement(current_coords, ccme);
      // element of current coordinate field w.r.t. reference coordinate mesh
      // elements
      ccrce = apf::createElement(current_coords, me);
    }
    void outElement()
    {
      apf::destroyElement(cccce);
      apf::destroyElement(ccrce);
      apf::destroyMeshElement(ccme);
      ElementalSystem::outElement();
    }
    bool includesBodyForces() { return true; }
    void atPoint(apf::Vector3 const &p, double w, double)
    {
      int &nen = nenodes;   // = 4 (tets)
      int &nedof = nedofs;  // = 12 (tets)
      apf::Matrix3x3 J, Jinv;
      // note that the jacobian is w.r.t the current coordinates!
      apf::getJacobian(ccme, p, J);
      apf::getJacobianInv(ccme, p, Jinv);
      double detJ = apf::getJacobianDeterminant(J, dim);
      apf::Matrix3x3 F;
      apf::getVectorGrad(ccrce, p, F);  // F=dx/dX
      double detF = getDeterminant(F);
      if (detF < 0.0) {
        std::cout << "error: detF < 0" << std::endl;
        exit(EXIT_FAILURE);
      }
      apf::NewArray<apf::Vector3> grads;
      apf::getShapeGrads(cccce, p, grads);
      // Calculate derivative of shape fxns wrt current coordinates.
      /*=======================================================
        Constitutive Component of Tangent Stiffness Matrix: K0
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        =======================================================
      */
      // hard-coded for 3d, make a general function... to produce this
      apf::DynamicMatrix BL(6, nedof);  // linear strain disp
      BL.zero();
      for (int ii = 0; ii < nen; ii++) {
        BL(0, dim * ii) = grads[ii][0];      // N_(ii,1)
        BL(1, dim * ii + 1) = grads[ii][1];  // N_(ii,2)
        BL(2, dim * ii + 2) = grads[ii][2];  // N_(ii,3)
        BL(3, dim * ii) = grads[ii][1];      // N_(ii,2)
        BL(3, dim * ii + 1) = grads[ii][0];  // N_(ii,1)
        BL(4, dim * ii + 1) = grads[ii][2];  // N_(ii,3)
        BL(4, dim * ii + 2) = grads[ii][1];  // N_(ii,2)
        BL(5, dim * ii) = grads[ii][2];      // N_(ii,3)
        BL(5, dim * ii + 2) = grads[ii][0];  // N_(ii,1)
      }
      /*
        Determine stress-strain matrix, D, for NeoHookean material.
        See Bonet and Wood 2nd Edition, PP.250.
        - lambda and mu are effective lame parameters.
        - lambda = ShearModulus/detJ, mu = ( PoissonsRatio - ShearModulus *
        ln(detJ) )/detJ
      */
      // Calculate Cauchy stress from compressible NeoHookean equation. See
      // Bonet and Wood 2nd Ed. PP.163.
      apf::setMatrix(dfm_grd_fld, apf::getMeshEntity(me),
                     current_integration_point, F);
      apf::DynamicMatrix leftCauchyGreen(3, 3);   // leftCauchyGreen.zero();
      apf::DynamicMatrix rightCauchyGreen(3, 3);  // rightCauchyGreen.zero();
      apf::DynamicMatrix FT(3, 3);
      FT.zero();
      apf::transpose(fromMatrix(F), FT);
      apf::multiply(fromMatrix(F), FT, leftCauchyGreen);
      apf::multiply(FT, fromMatrix(F), rightCauchyGreen);
      // apf::DynamicMatrix Cauchy(3,3);
      apf::Matrix3x3 Cauchy;
      double mu = ShearModulus;
      double lambda = (2.0 * mu * PoissonsRatio) / (1.0 - 2.0 * PoissonsRatio);
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          Cauchy[i][j] = mu / detF * (leftCauchyGreen(i, j) - (i == j)) +
                         lambda / detF * log(detF) * (i == j);
        }
      }
      double lambda_prime = lambda / detF;
      double mu_prime = (mu - lambda * log(detF)) / detF;
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
      apf::DynamicMatrix K0(nedof, nedof);
      apf::DynamicMatrix BLT(nedof, 6);
      apf::DynamicMatrix DBL(6, nedof);
      apf::transpose(BL, BLT);
      apf::multiply(D, BL, DBL);
      apf::multiply(BLT, DBL, K0);
      /*=======================================================================
        Initial Stress (or Geometric) component of Tangent Stiffness Matrix: K1
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        ========================================================================*/
      apf::DynamicMatrix BNL(9, nedof);  // nonlinear strain disp
      BNL.zero();
      for (int i = 0; i < nenodes; ++i) {
        BNL(0, i * dim) = grads[i][0];
        BNL(1, i * dim) = grads[i][1];
        BNL(2, i * dim) = grads[i][2];
        BNL(3, i * dim + 1) = grads[i][0];
        BNL(4, i * dim + 1) = grads[i][1];
        BNL(5, i * dim + 1) = grads[i][2];
        BNL(6, i * dim + 2) = grads[i][0];
        BNL(7, i * dim + 2) = grads[i][1];
        BNL(8, i * dim + 2) = grads[i][2];
      }
      apf::DynamicMatrix tau(9, 9);
      tau.zero();
      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
          tau(i, j) = Cauchy[i][j];
          tau(i + dim, j + dim) = Cauchy[i][j];
          tau(i + 2 * dim, j + 2 * dim) = Cauchy[i][j];
        }
      }
      apf::DynamicMatrix K1(nedof, nedof);
      K1.zero();
      apf::DynamicMatrix BNLT(nedof, 9);
      BNLT.zero();
      apf::DynamicMatrix BNLTxtau(9, nedof);
      BNLTxtau.zero();
      apf::transpose(BNL, BNLT);
      apf::multiply(BNLT, tau, BNLTxtau);
      apf::multiply(BNLTxtau, BNL, K1);
      /*=======================================================================
        Terms for force vector (RHS)
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        ========================================================================*/
      apf::DynamicVector cauchyVoigt(6);
      cauchyVoigt(0) = Cauchy[0][0];
      cauchyVoigt(1) = Cauchy[1][1];
      cauchyVoigt(2) = Cauchy[2][2];
      cauchyVoigt(3) = Cauchy[0][1];
      cauchyVoigt(4) = Cauchy[1][2];
      cauchyVoigt(5) = Cauchy[0][2];
      apf::DynamicVector BLTxCauchyVoigt(nedof);
      apf::multiply(BLT, cauchyVoigt, BLTxCauchyVoigt);
      for (int ii = 0; ii < nedof; ii++) {
        fe[ii] -= w * detJ * BLTxCauchyVoigt(ii);
        for (int jj = 0; jj < nedof; jj++)
          Ke(ii, jj) += w * detJ * (K0(ii, jj) + K1(ii, jj));
      }
      // E_G = 1/2(C-I), C=F^T.F
      apf::Matrix3x3 greenStrain(
          0.5 * (rightCauchyGreen(0, 0) - 1), 0.5 * rightCauchyGreen(0, 1),
          0.5 * rightCauchyGreen(0, 2), 0.5 * rightCauchyGreen(1, 0),
          0.5 * (rightCauchyGreen(1, 1) - 1), 0.5 * rightCauchyGreen(1, 2),
          0.5 * rightCauchyGreen(2, 0), 0.5 * rightCauchyGreen(2, 1),
          0.5 * (rightCauchyGreen(2, 2) - 1));
      analysis->storeStrain(me, greenStrain);
      analysis->storeStress(me, Cauchy);
      current_integration_point++;
    }
    int current_integration_point;

    private:
    NonlinearTissue *analysis;
    int dim;
    apf::Mesh *msh;
    apf::Field *dfm_grd_fld;
    // new stuff to try out new apf functions
    apf::Field *current_coords;
    apf::Element *cccce;
    apf::Element *ccrce;
    apf::MeshElement *ccme;
    //
    double ShearModulus;
    double PoissonsRatio;
  };
}
#endif
