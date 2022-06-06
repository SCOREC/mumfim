#ifndef MUMFIM_ULMULTISCALE_INTEGRATOR_H_
#define MUMFIM_ULMULTISCALE_INTEGRATOR_H_
#include <ElementalSystem.h>
#include <amsiDeformation.h>
#include <apfMatrixUtil.h>
#include <apfShape.h>
#include <cstring>
#include "NonlinearTissue.h"
#include "mumfim/macroscale/RVECoupling.h"
namespace mumfim
{
  // In this integrator we have made the assumption that the number of field
  // components is the dimenstion. We may want to replace the dim calls with
  // num_fild_components
  class ULMultiscaleIntegrator : public amsi::ElementalSystem
  {
    public:
    ULMultiscaleIntegrator(RVECoupling * r,
                           apf::Field * strn,
                           apf::Field * strs,
                           apf::Field * u,
                           apf::Field * dfm_grd,
                           int o,
                           double minimum_stiffness)
        : ElementalSystem(u, o)
        , current_integration_point(0)
        , coupling(r)
        , strain_field(strn)
        , stress_field(strs)
        , dfm_grd_fld(dfm_grd)
        , minimum_stiffness(minimum_stiffness)
    {
    }
    void inElement(apf::MeshElement * me) override
    {
      ElementalSystem::inElement(me);
      ref_lmnt =
          apf::createMeshElement(apf::getMesh(f), apf::getMeshEntity(me));
      du_lmnt = apf::createElement(f, ref_lmnt);
      dim = apf::getDimension(me);
    }
    void outElement() override
    {
      apf::destroyElement(du_lmnt);
      apf::destroyMeshElement(ref_lmnt);
      current_integration_point = 0;
      ElementalSystem::outElement();
    }
    bool includesBodyForces() override { return true; }
    // after the refactor, the micro-scale and macro-scale stress tesnors are
    // indexed identically Note that the system is processed using the current
    // coordinate field. I.E. the volume, and gradients of the default element
    // are updated to the current iteration
    void atPoint(apf::Vector3 const & p, double w, double dV) override
    {
      apf::MeshEntity * m = apf::getMeshEntity(me);
      micro_fo_result * rslt =
          coupling->getRVEResult(m, current_integration_point);
      int & nen = nenodes;   // = 4 (tets)
      int & nedof = nedofs;  // = 12 (tets)
      // fill the matrix form of the stress tensor
      apf::Matrix3x3 Cauchy;
      // This doesn't populate the matrix fully! only the upper diagonal!!!!
      // amsi::voigtVec2Mat(dim, &rslt->data[0], Cauchy);
      Cauchy[0][0] = rslt->data[0];
      Cauchy[1][1] = rslt->data[1];
      Cauchy[2][2] = rslt->data[2];
      Cauchy[1][2] = Cauchy[2][1] = rslt->data[3];
      Cauchy[0][2] = Cauchy[2][0] = rslt->data[4];
      Cauchy[0][1] = Cauchy[1][0] = rslt->data[5];
      apf::NewArray<apf::Vector3> grads;
      apf::getShapeGrads(e, p, grads);
      apf::NewArray<double> N;
      apf::getShapeValues(e, p, N);
      apf::DynamicMatrix BL(6, nedof);  // linear strain disp
      BL.zero();
      for (int ii = 0; ii < nen; ii++)
      {
        BL(0, dim * ii) = grads[ii][0];      // N_(ii,1)
        BL(1, dim * ii + 1) = grads[ii][1];  // N_(ii,2)
        BL(2, dim * ii + 2) = grads[ii][2];  // N_(ii,3)
        BL(5, dim * ii) = grads[ii][1];      // N_(ii,2)
        BL(5, dim * ii + 1) = grads[ii][0];  // N_(ii,1)
        BL(3, dim * ii + 1) = grads[ii][2];  // N_(ii,3)
        BL(3, dim * ii + 2) = grads[ii][1];  // N_(ii,2)
        BL(4, dim * ii) = grads[ii][2];      // N_(ii,3)
        BL(4, dim * ii + 2) = grads[ii][0];  // N_(ii,1)
      }
      apf::DynamicMatrix C(6, 6);
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
          C(i, j) = rslt->data[9 + i * 6 + j];
      
      apf::DynamicMatrix K0(nedof, nedof);
      apf::DynamicMatrix tmp(6, nedof);
      apf::DynamicMatrix BLT(nedof, 6);
      apf::transpose(BL, BLT);
      apf::multiply(C, BL, tmp);
      apf::multiply(BLT, tmp, K0);
      // Nonlinear terms (geometric stiffness)
      /*=======================================================================
        Initial Stress (or Geometric) component of Tangent Stiffness Matrix: K1
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        ========================================================================*/
      apf::DynamicMatrix BNL(9, nedof);  // nonlinear strain disp
      BNL.zero();
      for (int i = 0; i < nenodes; ++i)
      {
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
      for (int i = 0; i < dim; ++i)
      {
        for (int j = 0; j < dim; ++j)
        {
          tau(i, j) = Cauchy[i][j];
          tau(i + dim, j + dim) = Cauchy[i][j];
          tau(i + 2 * dim, j + 2 * dim) = Cauchy[i][j];
        }
      }
      apf::DynamicMatrix K1(nedof, nedof);
      K1.zero();
      apf::DynamicMatrix BNLT(nedof, 9);
      apf::DynamicMatrix BNLTxtau(9, nedof);
      apf::transpose(BNL, BNLT);
      apf::multiply(BNLT, tau, BNLTxtau);
      apf::multiply(BNLTxtau, BNL, K1);
      // for force vector calculation
      apf::DynamicMatrix BLTxSV(nedof, 1);
      BLTxSV.zero();
      // multiply linear portion of B with the Cauchy stress
      for (int ii = 0; ii < nedof; ii++)
        for (int jj = 0; jj < 6; jj++)
          BLTxSV(ii, 0) += BL(jj, ii) * rslt->data[jj];
      double wxdV = w * dV;
      for (int ii = 0; ii < nedof; ii++)
      {
        fe[ii] -= BLTxSV(ii, 0) * wxdV;  // P - F
        for (int jj = 0; jj < nedof; jj++)
        {
          Ke(ii, jj) += (K0(ii, jj) + K1(ii, jj)) * wxdV;
        }
      }
      // Compute the deformation gradient and strain tensor
      apf::Matrix3x3 F;
      amsi::deformationGradient(du_lmnt, p, F);
      apf::setMatrix(dfm_grd_fld, m, current_integration_point, F);
      // Calculate rightCauchyGreen tensor.
      apf::DynamicMatrix rightCauchyGreen(3, 3);  // rightCauchyGreen.zero();
      apf::DynamicMatrix FT(3, 3);
      FT.zero();
      apf::transpose(fromMatrix(F), FT);
      apf::multiply(FT, fromMatrix(F), rightCauchyGreen);
      // E_G = 1/2(C-I), C=F^T.F, Green-Lagrange Strain.
      apf::Matrix3x3 greenStrain(
          0.5 * (rightCauchyGreen(0, 0) - 1), 0.5 * rightCauchyGreen(0, 1),
          0.5 * rightCauchyGreen(0, 2), 0.5 * rightCauchyGreen(1, 0),
          0.5 * (rightCauchyGreen(1, 1) - 1), 0.5 * rightCauchyGreen(1, 2),
          0.5 * rightCauchyGreen(2, 0), 0.5 * rightCauchyGreen(2, 1),
          0.5 * (rightCauchyGreen(2, 2) - 1));
      apf::setMatrix(strain_field, m, current_integration_point, greenStrain);
      apf::setMatrix(stress_field, m, current_integration_point, Cauchy);
      current_integration_point++;
    }
    int current_integration_point;

    private:
    RVECoupling * coupling;
    apf::Field * strain_field;
    apf::Field * stress_field;
    int dim;
    apf::MeshElement * ref_lmnt;
    apf::Element * du_lmnt;
    apf::Field * dfm_grd_fld;
    double minimum_stiffness;
  };
}  // namespace mumfim
#endif
